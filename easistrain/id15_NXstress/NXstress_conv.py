import os
import h5py
import datetime
import tkinter as tk
from tkinter import ttk
import tkinter.font as tkFont
import numpy as np
import yaml
from easistrain.EDD.math import compute_qs
from easistrain.EDD.utils import fit_detector_data
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


def dt_str():
    return h5py.string_dtype(encoding="utf-8")


def _get_dtype(data):
    if isinstance(data, list) and all(isinstance(item, np.ndarray) for item in data):
        return h5py.special_dtype(vlen=np.float64)
    elif isinstance(data, np.ndarray):
        return data.dtype
    else:
        raise ValueError(
            "Data format not recognized. Expecting a numpy array or a list of numpy arrays."
        )


def _create_datasets(group, name, data):
    correct_dtype = _get_dtype(data)
    if len(data) == 1:
        dataset = group.create_dataset(name, data=data)
    else:
        dataset = group.create_dataset(name, data=data, dtype=correct_dtype)
    return dataset


def _create_nxtransformation(
    nxtransformations, name, value, units, ttype, vector, dtype=None
):
    if dtype is not None:
        value = np.array(value, dtype=dtype)
    else:
        value = np.array(value)

    dset = nxtransformations.create_dataset(name, data=value)
    if ttype is None:
        dset.attrs.update({"vector": vector})
    else:
        dset.attrs.update(
            {"transformation_type": ttype, "units": units, "vector": vector}
        )
    return dset


def _add_coordinate_system(nxtransformations, dset):
    dset.attrs["depends_on"] = "beam"
    dset = _create_nxtransformation(
        nxtransformations, "beam", float("nan"), None, None, [1, 0, 0]
    )
    dset.attrs["depends_on"] = "gravity"
    dset = _create_nxtransformation(
        nxtransformations, "gravity", float("nan"), None, None, [0, 0, -1]
    )
    dset.attrs["depends_on"] = "."


def _position_point_detector(nxdetector, distance, polar, azimuth):
    # X_lab =  Rx(-azimuth) . Ry(-polar) . Tx(distance) . X_sdd
    nxtransformations = nxdetector.create_group("position")
    nxtransformations.attrs["NX_class"] = "NXtransformations"
    nxdetector["depends_on"] = "position/distance"
    dset = _create_nxtransformation(
        nxtransformations, "distance", distance, "cm", "translation", [1, 0, 0]
    )
    dset.attrs["depends_on"] = "polar"
    dset = _create_nxtransformation(
        nxtransformations, "polar", polar, "deg", "rotation", [0, -1, 0]
    )
    dset.attrs["depends_on"] = "azimuth"
    dset = _create_nxtransformation(
        nxtransformations, "azimuth", azimuth, "deg", "rotation", [1, 0, 0]
    )
    _add_coordinate_system(nxtransformations, dset)


def _position_point_detector_slits(nxslit, distance, polar, azimuth):
    # X_lab =  Rx(azimuth) . Ry(-polar) . Tx(distance) . Ry(90) . X_sdd
    nxtransformations = nxslit.create_group("position")
    nxtransformations.attrs["NX_class"] = "NXtransformations"
    nxslit["depends_on"] = "position/perpendicular"
    dset = _create_nxtransformation(
        nxtransformations,
        "perpendicular",
        90,
        "deg",
        "rotation",
        [0, 1, 0],
        dtype=np.int64,
    )
    dset.attrs["depends_on"] = "distance"
    dset = _create_nxtransformation(
        nxtransformations, "distance", distance, "cm", "translation", [1, 0, 0]
    )
    dset.attrs["depends_on"] = "polar"
    dset = _create_nxtransformation(
        nxtransformations, "polar", polar, "deg", "rotation", [0, -1, 0]
    )
    dset.attrs["depends_on"] = "azimuth"
    dset = _create_nxtransformation(
        nxtransformations, "azimuth", azimuth, "deg", "rotation", [1, 0, 0]
    )
    _add_coordinate_system(nxtransformations, dset)


def _position_source_slits(nxslit):
    # X_lab = Tx(-distance) . Ry(90) . X_sdd
    nxtransformations = nxslit.create_group("position")
    nxtransformations.attrs["NX_class"] = "NXtransformations"

    nxslit["depends_on"] = "position/perpendicular"
    dset = _create_nxtransformation(
        nxtransformations,
        "perpendicular",
        90,
        "deg",
        "rotation",
        [0, 1, 0],
        dtype=np.int64,
    )
    dset.attrs["depends_on"] = "distance"
    dset = _create_nxtransformation(
        nxtransformations, "distance", float("nan"), "cm", "translation", [-1, 0, 0]
    )
    _add_coordinate_system(nxtransformations, dset)


class NXstressFromRaw:
    def __init__(
        self,
        file_path,
        det_calib_file_energy,
        det_calib_file_angle,
        lattice,
        phase_name,
        experimental_identifier,
        collection_identifier,
        with_cradle=False,
        scanNbForRotation=None,
        test_script=False,
    ):

        self.file_path = file_path
        self.energy_calib_det0, self.energy_calib_det1 = self.extract_energy_calib(
            det_calib_file_energy
        )
        self.energy_calib = None

        self.lattice = lattice
        self.phase_name = phase_name
        self.experimental_identifier = experimental_identifier
        self.collection_identifier = collection_identifier

        self.diffractogram = None
        self.with_cradle = with_cradle
        self.scanNbForRotation = scanNbForRotation  # The index of the scan in raw file where sample was rotated
        self.rotated = False

        self.rangeFit = None
        self.rangeFitDet1 = None
        self.rangeFitDet2 = (
            None  # This really just detector 0, but for when we rotate the sample
        )
        self.nbPeaksInBoxes = []
        self.nbPeaksInBoxesDet1 = []
        self.nbPeaksInBoxesDet2 = []  # This is the same for self.rangeFitDet2

        self.results = {}
        self.diffractogram_count = 0

        self.phi = None
        self.chi = None
        self.omega = None

        self.mca_det0_polar, self.mca_det1_polar = self.extract_angle_calib(
            det_calib_file_angle
        )

        self.mca_det0_azimuth = -90.0  # I'm guessing this is always the case?
        self.mca_det1_azimuth = 0.0

        self.mca_det0_distance = float("nan")
        self.mca_det1_distance = float("nan")

        self.h_miller = []
        self.k_miller = []
        self.l_miller = []

        if test_script is True:
            self.h_miller = [3]
            self.k_miller = [1]
            self.l_miller = [1]
            self.rangeFit = [(1651, 1949)]
            self.nbPeaksInBoxes = [2]
            self.rangeFitDet1 = [(1601, 1871)]
            self.nbPeaksInBoxesDet1 = [2]
            self.rangeFitDet2 = [(1659, 1958)]
            self.nbPeaksInBoxesDet2 = [2]

    def extract_angle_calib(self, file_path):
        with h5py.File(file_path, "r") as file:
            mca_det0_polar = self.find_dataset(file, "calibratedAngleHD")
            mca_det1_polar = self.find_dataset(file, "calibratedAngleVD")
            return np.float64(mca_det0_polar), np.float64(mca_det1_polar)

    def extract_energy_calib(self, file_path):
        with h5py.File(file_path, "r") as file:
            mca_det0_calib = self.find_dataset(file, "calibCoeffsHD")
            mca_det1_calib = self.find_dataset(file, "calibCoeffsVD")
            return mca_det0_calib, mca_det1_calib

    def find_dataset(self, group, dataset_name):
        """
        Recursively searches all subgroups of the given group for the dataset.
        Returns the dataset if found, or None if not found.
        """
        if dataset_name in group:
            return group[dataset_name][:]
        for sub_group_name in group:
            sub_group = group[sub_group_name]
            if isinstance(sub_group, h5py.Group):
                found_dataset = self.find_dataset(sub_group, dataset_name)
                if found_dataset is not None:
                    return found_dataset
        return None

    def create_beam_intensity_profile(self, output_file, entry_name):
        """
        Creates a 'beam_intensity_profile' subgroup within a specified entry in an HDF5 file.
        This subgroup is designed to store detailed settings and parameters that define the beam
        intensity profile used during the experiment.

        Parameters:
        - output_file: The HDF5 file object open for writing.
        - entry_name: The name of the entry under which the beam intensity profile group is created.

        The method organizes the subgroup according to the NeXus NXcollection class, documenting
        various parameters related to primary and secondary beam shaping elements such as slits.
        Attributes such as type, full width, width, and distance for both primary and secondary
        vertical and horizontal beam components are defined, ensuring precise documentation of
        the experimental setup. This structured data supports reproducibility and detailed analysis
        of experimental conditions.
        """
        entry_group = output_file.require_group(entry_name)
        beam_intensity_profile_group = entry_group.create_group(
            "beam_intensity_profile"
        )

        beam_intensity_profile_group.attrs["NX_class"] = "NXcollection"

        # Adding gauge volume
        parameters = {
            "primary_vertical_type": "slit",
            "primary_vertical_full_width": 100.0,
            "primary_vertical_width": 100.0,
            "primary_vertical_distance": float("nan"),
            "primary_horizontal_type": "slit",
            "primary_horizontal_full_width": 100.0,
            "primary_horizontal_width": 100.0,
            "primary_horizontal_distance": float("nan"),
            "secondary_horizontal_type": "slit",
            "secondary_horizontal_full_width": 100.0,
            "secondary_horizontal_width": 100.0,
            "secondary_horizontal_distance": float("nan"),
        }

        for parameter_name, parameter_value in parameters.items():
            dataset = beam_intensity_profile_group.create_dataset(
                parameter_name, data=parameter_value
            )
            if "width" in parameter_name:
                dataset.attrs["units"] = "um"
            elif "distance" in parameter_name:
                dataset.attrs["units"] = "mm"

    def create_fit_subgroup(self, output_file, entry_name):
        """
        Creates a 'fit' subgroup within a specified entry in an HDF5 file, organizing and storing various
        fit-related data and parameters processed during a diffraction analysis. This subgroup includes
        detailed datasets for fit parameters, background parameters, and a comprehensive description
        of the fitting process.

        Parameters:
        - output_file: The HDF5 file object open for writing.
        - entry_name: The name of the entry under which the fit group is created.

        The method structures the subgroup according to the NeXus NXprocess class, defining data for
        each fit parameter along with its uncertainty. It captures the diffractogram, associated fitted
        data, and background data, along with residual calculations. It also manages data on peak
        analysis parameters, such as area, center, and full width at half maximum (FWHM), ensuring
        that all parameters and their errors are accurately documented. The organization and labeling
        of the data adhere to standards that facilitate subsequent data retrieval and analysis.
        """

        def process_fit_parameters(data):
            columns = [
                "area",
                "center",
                "fwhm_left",
                "fwhm_right",
                "form_factor",
                "goodness_of_fit",
            ]
            data = np.asarray(data).T
            fitParams = dict(zip(columns, data))
            return fitParams

        entry_group = output_file[entry_name]
        fit_group = entry_group.create_group("fit")

        fit_group.attrs["NX_class"] = "NXprocess"
        fit_group.attrs["default"] = "diffractogram"
        fit_group.create_dataset("version", data="1.0")
        fit_group.create_dataset(
            "date", data=datetime.datetime.now().strftime("%Y-%m-%d")
        )
        # fit_group.create_dataset('sequence_index', data=[1], dtype=np.int64)
        fit_group.create_dataset("program", data="easistrain")

        background_parameters_group = fit_group.create_group("background_parameters")
        background_parameters_group.attrs["NX_class"] = "NXdata"
        background_parameters_group.attrs["auxiliary_signals"] = ["A1"]
        background_parameters_group.attrs["signal"] = "A0"
        background_parameters_group.create_dataset("A0", data=np.NaN)
        background_parameters_group.create_dataset("A1", data=np.NaN)
        background_parameters_group.create_dataset("title", data="linear")

        description_group = fit_group.create_group("description")
        description_group.attrs["NX_class"] = "NXnote"
        description_group.attrs["type"] = "text/plain"
        description_group.create_dataset("data", data="fitted by easistrain 1.0")

        diffractogram_group = fit_group.create_group("diffractogram")
        diffractogram_group.attrs["NX_class"] = "NXdata"
        diffractogram_group.attrs["auxiliary_signals"] = [
            "fit",
            "background",
            "residuals",
        ]
        diffractogram_group.attrs["axes"] = ["patterns", "channels"]
        diffractogram_group.attrs["signal"] = "diffractogram"
        diffractogram_group.create_dataset("diffractogram", data=self.diffractogram)
        _create_datasets(diffractogram_group, "background", self.results["background"])
        _create_datasets(diffractogram_group, "fit", self.results["fittedData"])

        channels = []
        for i, (start, end) in enumerate(self.rangeFit):
            channels.append(np.arange(start, end))
        _create_datasets(diffractogram_group, "channels", channels)

        peak_parameters_group = fit_group.create_group("peak_parameters")
        peak_parameters_group.attrs["NX_class"] = "NXdata"
        peak_parameters_group.attrs["signal"] = "area"
        peak_parameters_group.attrs["auxiliary_signals"] = [
            "center",
            "fwhm_left",
            "fwhm_right",
            "form_factor",
            "goodness_of_fit",
        ]

        fitParams = process_fit_parameters(self.results["fitParams"])
        area_dataset = peak_parameters_group.create_dataset(
            "area", data=fitParams["area"]
        )
        area_dataset.attrs["units"] = "channels"
        center_dataset = peak_parameters_group.create_dataset(
            "center", data=fitParams["center"]
        )
        center_dataset.attrs["units"] = "channels"
        fwhm_left_dataset = peak_parameters_group.create_dataset(
            "fwhm_left", data=fitParams["fwhm_left"]
        )
        fwhm_left_dataset.attrs["units"] = "channels"
        fwhm_right_dataset = peak_parameters_group.create_dataset(
            "fwhm_right", data=fitParams["fwhm_right"]
        )
        fwhm_right_dataset.attrs["units"] = "channels"
        peak_parameters_group.create_dataset(
            "form_factor", data=fitParams["form_factor"]
        )
        peak_parameters_group.create_dataset(
            "goodness_of_fit", data=fitParams["goodness_of_fit"]
        )

        uncertaintyfitParams = process_fit_parameters(
            self.results["uncertaintyFitParams"]
        )
        uncertainty_area_dataset = peak_parameters_group.create_dataset(
            "area_error", data=uncertaintyfitParams["area"]
        )
        uncertainty_area_dataset.attrs["units"] = "channels"
        uncertainty_center_dataset = peak_parameters_group.create_dataset(
            "center_error", data=uncertaintyfitParams["center"]
        )
        uncertainty_center_dataset.attrs["units"] = "channels"
        uncertainty_fwhm_left_dataset = peak_parameters_group.create_dataset(
            "fwhm_left_error", data=uncertaintyfitParams["fwhm_left"]
        )
        uncertainty_fwhm_left_dataset.attrs["units"] = "channels"
        uncertainty_fwhm_right_dataset = peak_parameters_group.create_dataset(
            "fwhm_right_error", data=uncertaintyfitParams["fwhm_right"]
        )
        uncertainty_fwhm_right_dataset.attrs["units"] = "channels"
        peak_parameters_group.create_dataset(
            "form_factor_error", data=uncertaintyfitParams["form_factor"]
        )

    def create_instrument_subgroup(self, output_file, entry_name):
        """
        Constructs an instrument subgroup in the HDF5 output file under a specified entry.
        This subgroup contains detailed descriptions and settings for the instrument components,
        including detectors, slits, and the source, all formatted according to NeXus standards.

        Parameters:
        - output_file: The HDF5 file object open for writing.
        - entry_name: The name of the entry under which the instrument group is created.

        This method sets up the structure and populates it with information regarding the instrument
        configuration used during the experiment, ensuring that all components like detectors,
        slits, and source are accurately represented.
        """
        entry_group = output_file[entry_name]
        instrument_group = entry_group.create_group("instrument")

        instrument_group.attrs["NX_class"] = "NXinstrument"
        instrument_group.create_dataset("name", data="ESRF-ID15")

        mca2_det0_group = instrument_group.create_group("mca2_det0")
        mca2_det0_group.attrs["NX_class"] = "NXdetector"
        mca2_det0_group["type"] = "MCA"
        mca2_det0_group["description"] = "Ge MCA"
        _position_point_detector(
            mca2_det0_group,
            self.mca_det0_distance,
            self.mca_det0_polar,
            self.mca_det0_azimuth,
        )

        mca2_det0_slits_group = instrument_group.create_group("mca2_det0_slits")
        mca2_det0_slits_group.attrs["NX_class"] = "NXslit"
        mca2_det0_slits_group["x_gap"] = 100.0
        mca2_det0_slits_group["x_gap"].attrs["units"] = "um"
        mca2_det0_slits_group["y_gap"] = 100.0
        mca2_det0_slits_group["y_gap"].attrs["units"] = "um"
        _position_point_detector_slits(
            mca2_det0_slits_group,
            self.mca_det0_distance,
            self.mca_det0_polar,
            self.mca_det0_azimuth,
        )

        mca2_det1_group = instrument_group.create_group("mca2_det1")
        mca2_det1_group.attrs["NX_class"] = "NXdetector"
        mca2_det1_group["type"] = "MCA"
        mca2_det1_group["description"] = "Ge MCA"
        _position_point_detector(
            mca2_det1_group,
            self.mca_det1_distance,
            self.mca_det1_polar,
            self.mca_det1_azimuth,
        )

        mca2_det1_slits_group = instrument_group.create_group("mca2_det1_slits")
        mca2_det1_slits_group.attrs["NX_class"] = "NXslit"
        mca2_det1_slits_group["x_gap"] = 100.0
        mca2_det1_slits_group["x_gap"].attrs["units"] = "um"
        mca2_det1_slits_group["y_gap"] = 100.0
        mca2_det1_slits_group["y_gap"].attrs["units"] = "um"
        _position_point_detector_slits(
            mca2_det1_slits_group,
            self.mca_det1_distance,
            self.mca_det1_polar,
            self.mca_det1_azimuth,
        )

        primary_slits_group = instrument_group.create_group("primary_slits")
        primary_slits_group.attrs["NX_class"] = "NXslit"
        primary_slits_group["x_gap"] = 100.0
        primary_slits_group["x_gap"].attrs["units"] = "um"
        primary_slits_group["y_gap"] = 100.0
        primary_slits_group["y_gap"].attrs["units"] = "um"
        _position_source_slits(primary_slits_group)

        source_group = instrument_group.create_group("source")
        source_group.attrs["NX_class"] = "NXsource"
        source_group.create_dataset("probe", data="X-ray")
        source_group.create_dataset("type", data="Synchrotron X-ray Source")

    def create_peaks_subgroup(
        self, output_file, entry_name, sx_value, sy_value, sz_value, detector_index
    ):
        """
        Creates and populates a 'peaks' subgroup within a given entry in an HDF5 file.
        This subgroup stores detailed diffraction and position data for one or more peaks,
        formatted according to the NeXus NXdata standard. This method supports analyses
        that involve multiple peaks, dynamically adapting to the number of peaks specified.

        Parameters:
        - output_file: The HDF5 file object open for writing.
        - entry_name: The name of the entry under which the peaks group is created.
        - sx_value: The sx position value for the peak(s). Supports both single and multiple values.
        - sy_value: The sy position value for the peak(s). Supports both single and multiple values.
        - sz_value: The sz position value for the peak(s). Supports both single and multiple values.
        - detector_index: Index specifying which detector's settings to use for azimuth and polar angles.

        This method ensures that each peak's spatial and diffraction data are correctly recorded,
        including auxiliary signals and transformation metadata, while accommodating setups with
        multiple peak analyses.
        """
        if entry_name not in output_file:
            entry_group = output_file.create_group(entry_name)
        else:
            entry_group = output_file[entry_name]

        peaks_group = entry_group.require_group("peaks")
        peaks_group.attrs["NX_class"] = "NXdata"

        peaks_group.attrs["signal"] = "h"
        peaks_group.attrs["auxiliary_signals"] = [
            "k",
            "l",
            "qx",
            "qy",
            "qz",
            "sx",
            "sy",
            "sz",
            "center",
        ]

        if detector_index == 0:
            det_azim = self.mca_det0_azimuth
            det_polar = self.mca_det0_polar
        elif detector_index == 1:
            det_azim = self.mca_det1_azimuth
            det_polar = self.mca_det1_polar
        else:
            raise ValueError(
                "Invalid detector index. Only 0 (horizontal) and 1 (vertical) are valid."
            )

        angles = [self.phi, self.chi, self.omega, det_azim, det_polar]

        n = len(self.nbPeaksInBoxes)

        qx, qy, qz = compute_qs(angles)

        if self.rotated:
            qx, qy = qy, qx

        if n > 1:
            qx = np.full(n, qx)
            qy = np.full(n, qy)
            qz = np.full(n, qz)
            sx_value = np.full(n, sx_value)
            sy_value = np.full(n, sy_value)
            sz_value = np.full(n, sz_value)

        peaks_group.create_dataset("qx", data=qx, dtype=np.float64)
        peaks_group.create_dataset("qy", data=qy, dtype=np.float64)
        peaks_group.create_dataset("qz", data=qz, dtype=np.float64)

        peaks_group.create_dataset("sx", data=sx_value, dtype=np.float64)
        peaks_group.create_dataset("sy", data=sy_value, dtype=np.float64)
        peaks_group.create_dataset("sz", data=sz_value, dtype=np.float64)

        n = len(self.results["fitParams"])

        center = np.asarray(
            [self.results["fitParams"][i][1] for i in range(n)], dtype=np.float64
        )
        center_error = np.asarray(
            [self.results["uncertaintyFitParams"][i][1] for i in range(n)],
            dtype=np.float64,
        )

        energy = (
            center**2 * self.energy_calib[0]
            + center * self.energy_calib[1]
            + self.energy_calib[2]
        )

        energy_error = (
            np.sqrt(
                (2 * center * self.energy_calib[0]) ** 2 + (self.energy_calib[1]) ** 2
            )
            * center_error
        )

        peaks_group.create_dataset("center", data=energy)
        peaks_group["center"].attrs["units"] = "keV"
        peaks_group.create_dataset("center_errors", data=energy_error)
        peaks_group["center_errors"].attrs["units"] = "keV"

        h_data = np.asarray(self.h_miller, dtype=np.int64)
        k_data = np.asarray(self.k_miller, dtype=np.int64)
        l_data = np.asarray(self.l_miller, dtype=np.int64)

        peaks_group.create_dataset("h", data=h_data)
        peaks_group.create_dataset("k", data=k_data)
        peaks_group.create_dataset("l", data=l_data)

        peaks_group.create_dataset("center_type", data="energy", dtype=dt_str())
        peaks_group.create_dataset("lattice", data=self.lattice, dtype=dt_str())
        peaks_group.create_dataset("phase_name", data=self.phase_name, dtype=dt_str())
        peaks_group.create_dataset("title", data="peak parameters", dtype=dt_str())

    def perform_fit(self):
        """
        Performs fitting operations on selected ranges of the diffractogram. Each range is analyzed
        to extract peak parameters and background based on predefined settings and the number of peaks
        to be fitted within each segment.

        Returns:
        - results: A dictionary containing arrays of fit parameters, uncertainties, background data,
                and the actual fitted data for each processed range.

        This method leverages the `fit_detector_data` function to apply fitting procedures to the segments
        of the diffractogram specified by `rangeFit`.
        """
        results = {
            "fitParams": [],
            "uncertaintyFitParams": [],
            "background": [],
            "fittedData": [],
        }

        for i, (fit_min, fit_max) in enumerate(self.rangeFit):
            nb_peaks = self.nbPeaksInBoxes[i]
            channels = np.arange(fit_min, fit_max)
            raw_data = self.diffractogram[fit_min:fit_max]
            assert isinstance(raw_data, np.ndarray)

            (
                background,
                fitted_data,
                fit_params,
                uncertainty_fit_params,
            ) = fit_detector_data(
                channels=channels,
                raw_data=raw_data,
                nb_peaks=nb_peaks,
                boxCounter=i,
                scanNumber=0,
                detectorName="",
            )

            results["background"].append(background)
            results["fittedData"].append(fitted_data)
            results["fitParams"].append(fit_params)
            results["uncertaintyFitParams"].append(uncertainty_fit_params)

        return results

    def create_sample_subgroup(
        self, output_file, raw_file_entry, entry_name, sx_value, sy_value, sz_value
    ):
        """
        Creates a sample subgroup within the given entry in the HDF5 output file,
        defining the sample's position and orientation transformations.

        Args:
        - output_file: The HDF5 file handle where the subgroup is to be created.
        - raw_file_entry: The entry from the raw data file used as a source for sample information.
        - entry_name: The name of the entry under which the subgroup is created.
        - sx_value: The x-coordinate value of the sample's position.
        - sy_value: The y-coordinate value of the sample's position.
        - sz_value: The z-coordinate value of the sample's position.

        This method sets up a structured way of storing transformation data for the sample
        related to its position and orientation, following the NXsample and NXtransformations classes.
        """
        entry_group = output_file[entry_name]
        sample_group = entry_group.create_group("sample")

        sample_group.attrs["NX_class"] = "NXsample"

        nxtransformations = sample_group.create_group("position")
        nxtransformations.attrs["NX_class"] = "NXtransformations"

        sample_group["depends_on"] = "position/phi"
        dset = _create_nxtransformation(
            nxtransformations,
            "phi",
            self.phi,
            "deg",
            "rotation",
            [0, 0, -1],
        )
        dset.attrs["depends_on"] = "chi"
        dset = _create_nxtransformation(
            nxtransformations,
            "chi",
            self.chi,
            "deg",
            "rotation",
            [1, 0, 0],
        )
        dset.attrs["depends_on"] = "omega"
        dset = _create_nxtransformation(
            nxtransformations,
            "omega",
            self.omega,
            "deg",
            "rotation",
            [0, -1, 0],
        )
        dset.attrs["depends_on"] = "x"
        dset = _create_nxtransformation(
            nxtransformations,
            "x",
            sx_value,
            "mm",
            "translation",
            [1, 0, 0],
        )
        dset.attrs["depends_on"] = "y"
        dset = _create_nxtransformation(
            nxtransformations,
            "y",
            sy_value,
            "mm",
            "translation",
            [0, 1, 0],
        )
        dset.attrs["depends_on"] = "z"
        dset = _create_nxtransformation(
            nxtransformations,
            "z",
            sz_value,
            "mm",
            "translation",
            [0, 0, 1],
        )
        _add_coordinate_system(nxtransformations, dset)

        with h5py.File(self.file_path, "r") as raw_file:
            sample_name = raw_file[raw_file_entry]["start_time"]
            sample_group.create_dataset(
                "name", data=sample_name[()], dtype=sample_name.dtype
            )

    def create_entry_datasets(self, output_file, raw_file_entry, entry_name):
        """
        Creates and populates datasets within a new entry in the output HDF5 file.
        This method sets various attributes and datasets to comply with the NeXus NXstress standard for
        data representation, including time stamps, experiment identifiers, and data types.

        Args:
        - output_file: The HDF5 file object where data should be written.
        - raw_file_entry: The entry from the raw data file used as a source for some datasets.
        - entry_name: The name of the new entry being created in the output file.

        The method initializes an entry group and populates it with metadata and data relevant
        to the stress analysis experiments conducted.
        """
        entry_group = output_file[entry_name]

        entry_group.attrs["NX_class"] = "NXentry"
        entry_group.attrs["default"] = "fit"

        dt_str = h5py.string_dtype(encoding="utf-8")

        entry_group.create_dataset("definition", data="NXstress", dtype=dt_str)

        with h5py.File(self.file_path, "r") as raw_file:
            start_time = raw_file[raw_file_entry]["start_time"]
            end_time = raw_file[raw_file_entry]["end_time"]
            entry_group.create_dataset(
                "start_time", data=start_time[()], dtype=start_time.dtype
            )
            entry_group.create_dataset(
                "end_time", data=end_time[()], dtype=end_time.dtype
            )

        entry_group.create_dataset("diffraction_type", data="energy", dtype=dt_str)

        entry_group.create_dataset(
            "experiment_identifier", data=self.experimental_identifier, dtype=dt_str
        )
        entry_group.create_dataset(
            "experiment_description",
            data="energy-dispersive X-ray powder diffraction",
            dtype=dt_str,
        )
        entry_group.create_dataset(
            "collection_identifier", data=self.collection_identifier, dtype=dt_str
        )

    def initial_prompts(self, raw_file, raw_file_entry, idx):
        """
        Handles initial data prompting for peak and hkl value selection based on detector data.
        Launches GUI prompts to set peak fitting ranges and hkl values if needed.

        Args:
        - raw_file: The HDF5 file handle containing the measurement data.
        - raw_file_entry: The specific entry in the HDF5 file to be processed.
        - idx: The index representing the specific dataset within the entry to process.
        """

        def handle_gui_data(peak_ranges, peaks, detector_id):
            """
            Handles the GUI data submission for setting peak ranges.
            Updates class attributes based on the detector ID.

            Args:
            - peak_ranges: The selected range of peaks from the GUI.
            - peaks: String of comma-separated peak values.
            - detector_id: The identifier for the detector (0, 1, or 2).
            """
            if detector_id == 0:
                self.rangeFit = peak_ranges
                self.nbPeaksInBoxes.extend(map(int, peaks.split(",")))
                print(
                    f"Updated for detector 0 - rangeFit: {self.rangeFit}, nbPeaksInBoxes: {self.nbPeaksInBoxes}"
                )

            elif detector_id == 1:
                self.rangeFitDet1 = peak_ranges
                self.nbPeaksInBoxesDet1.extend(map(int, peaks.split(",")))
                print(
                    f"Updated for detector 1 - rangeFitDet1: {self.rangeFitDet1}, nbPeaksInBoxesDet1: {self.nbPeaksInBoxesDet1}"
                )

            elif detector_id == 2:
                self.rangeFitDet2 = peak_ranges
                self.nbPeaksInBoxesDet2.extend(map(int, peaks.split(",")))
                print(
                    f"Updated for detector 2 - rangeFitDet2: {self.rangeFitDet2}, nbPeaksInBoxesDet2: {self.nbPeaksInBoxesDet2}"
                )

        def handle_hkl_data(h_miller, k_miller, l_miller):
            """
            Handles the GUI data submission for hkl values.
            Updates class attributes to store hkl values from the GUI.

            Args:
            - h_miller: List of h-values.
            - k_miller: List of k-values.
            - l_miller: List of l-values.
            """
            self.h_miller.extend(h_miller)
            self.k_miller.extend(k_miller)
            self.l_miller.extend(l_miller)

        if self.rangeFit is None:
            detector0_path = f"{raw_file_entry}/measurement/mca2_det0"
            if idx < raw_file[detector0_path].shape[0]:
                diffractogram0 = raw_file[detector0_path][idx]
                if np.max(diffractogram0) > 1000:
                    gui = GUI(
                        diffractogram0,
                        lambda ranges, peaks: handle_gui_data(ranges, peaks, 0),
                        title="Select fit ranges for detector 0",
                    )
                    gui.mainloop()

        if self.rangeFitDet1 is None:
            start_index = 0
            while True:
                detector1_path = f"{raw_file_entry}/measurement/mca2_det1"
                if (
                    detector1_path in raw_file
                    and start_index < raw_file[detector1_path].shape[0]
                ):
                    diffractogram1 = raw_file[detector1_path][start_index]
                    if np.max(diffractogram1) > 1000:
                        gui = GUI(
                            diffractogram1,
                            lambda ranges, peaks: handle_gui_data(ranges, peaks, 1),
                            title="Select fit ranges for detector 1",
                        )
                        gui.mainloop()
                        break
                start_index += 1
                if start_index >= len(raw_file):
                    break

        if self.scanNbForRotation is not None and self.rangeFitDet2 is None:
            start_index = int(self.scanNbForRotation) + 1
            while True:
                entry = f"{start_index}.1"
                path = f"{entry}/measurement/mca2_det0"
                if (
                    path in raw_file
                    and idx < raw_file[path].shape[0]
                    and np.max(raw_file[path][idx]) > 1000
                ):
                    diffractogram2 = raw_file[path][idx]
                    gui = GUI(
                        diffractogram2,
                        lambda ranges, peaks: handle_gui_data(ranges, peaks, 2),
                        title="Select fit ranges for detector 2",
                    )
                    gui.mainloop()
                    break
                start_index += 1

        if (
            len(self.h_miller) == 0
            and len(self.k_miller) == 0
            and len(self.l_miller) == 0
        ):
            gui_hkl = GUI(
                data=[1],
                handle_gui_data=handle_hkl_data,
                title="Prompt hkl values",
                prompt_hkl=True,
            )
            gui_hkl.mainloop()

    def get_positioners(self, raw_file, raw_file_entry):
        """
        Extracts the positioner values from a specific entry in the provided HDF5 file.
        Assumes raw_file is an open HDF5 file.

        Args:
            raw_file: Opened HDF5 file object.
            raw_file_entry (str): The entry key within the HDF5 file.

        Returns:
            tuple: Contains the positioner values (sx, sy, sz, phi, chi & omega).
        """
        if self.with_cradle:
            phi = raw_file[raw_file_entry]["instrument/positioners/ephi"][()]
            chi = raw_file[raw_file_entry]["instrument/positioners/echi"][()]
            omega = raw_file[raw_file_entry]["instrument/positioners/sy"][()]
            sx = raw_file[raw_file_entry]["instrument/positioners/ex"][()]
            sy = raw_file[raw_file_entry]["instrument/positioners/ey"][()]
            sz = raw_file[raw_file_entry]["instrument/positioners/ez"][()]

        else:
            phi = raw_file[raw_file_entry]["instrument/positioners/srz"][()]
            chi = raw_file[raw_file_entry]["instrument/positioners/ssy2"][()]
            omega = raw_file[raw_file_entry]["instrument/positioners/y42"][()]
            sx = raw_file[raw_file_entry]["instrument/positioners/sx"][()]
            sy = raw_file[raw_file_entry]["instrument/positioners/sy"][()]
            sz = raw_file[raw_file_entry]["instrument/positioners/sz"][()]

        if self.rotated:
            sx, sy = sy, sx

        return sx, sy, sz, phi, chi, omega

    def check_detector(self, raw_file, raw_file_entry, idx, detector_index):
        """
        Checks if the detector data at a specified index passes a certain threshold criterion.
        This method is used to verify if the maximum value in a diffractogram is above a defined limit (1000),
        indicating that the data is significant enough for further processing.

        Args:
        - raw_file: The HDF5 file handle containing the measurement data.
        - raw_file_entry: The specific entry in the HDF5 file to be processed.
        - idx: The index within the entry indicating the specific dataset to check.
        - detector_index: The index of the detector being checked.

        Returns:
        - passes_check: A boolean value indicating whether the data passes the threshold check.
        """
        passes_check = True
        try:
            detector_path = f"{raw_file_entry}/measurement/mca2_det{detector_index}"
            if idx < raw_file[detector_path].shape[0]:
                diffractogram = raw_file[detector_path][idx]
                if np.max(diffractogram) <= 1000:
                    passes_check = False
        except KeyError:
            passes_check = False
            print(
                f"Missing dataset for detector {detector_index} in entry '{raw_file_entry}'. Skipping this entry."
            )

        return passes_check

    def process_scan(self, output_file_name, raw_file_entry, detector_index):
        """
        Processes the specified HDF5 file to evaluate and record diffractogram data for given positioners and detectors.
        This method also initiates data processing such as peak fitting and intensity profile creation if certain conditions are met.

        Args:
        - output_file_name: Name of the output HDF5 file where results are stored.
        - raw_file_entry: The entry in the HDF5 file corresponding to the dataset being processed.
        - detector_index: Index of the detector used to fetch the diffractogram data.

        The method checks if the diffractogram passes a specified value check, prompts for additional input if necessary,
        and creates structured entries in an output file based on the processed data.
        """
        with h5py.File(self.file_path, "r") as raw_file:
            sx, sy, sz, self.phi, self.chi, self.omega = self.get_positioners(
                raw_file, raw_file_entry
            )

            positioners = ["sx", "sy", "sz"]
            _, pos_array = None, None
            for pos in positioners:
                data_path = f"instrument/positioners/{pos}"
                if data_path in raw_file[raw_file_entry]:
                    data = raw_file[raw_file_entry][data_path][()]
                    if isinstance(data, np.ndarray) and data.size > 1:
                        _, pos_array = pos, data
                        break

            n_positions = len(pos_array) if pos_array is not None else 1

            with h5py.File(output_file_name, "a") as output_file:
                for idx in range(n_positions):
                    sx_val = sx[idx] if isinstance(sx, np.ndarray) else sx
                    sy_val = sy[idx] if isinstance(sy, np.ndarray) else sy
                    sz_val = sz[idx] if isinstance(sz, np.ndarray) else sz

                    passes_check = self.check_detector(
                        raw_file, raw_file_entry, idx, detector_index
                    )

                    if passes_check:
                        try:
                            if self.rangeFit is None:
                                self.initial_prompts(raw_file, raw_file_entry, idx)
                            detector_path = (
                                f"{raw_file_entry}/measurement/mca2_det{detector_index}"
                            )
                            if idx < raw_file[detector_path].shape[0]:
                                self.diffractogram = raw_file[detector_path][idx]
                        except KeyError:
                            print(
                                f"Missing dataset for detector {detector_index} in entry '{raw_file_entry}'. Skipping prompt for rangeFit."
                            )

                        self.diffractogram_count += 1
                        modified_entry_name = f"{self.diffractogram_count}.1"

                        self.results = self.perform_fit()

                        self.create_peaks_subgroup(
                            output_file,
                            modified_entry_name,
                            sx_val,
                            sy_val,
                            sz_val,
                            detector_index,
                        )
                        self.create_beam_intensity_profile(
                            output_file, modified_entry_name
                        )
                        self.create_fit_subgroup(output_file, modified_entry_name)
                        self.create_instrument_subgroup(
                            output_file, modified_entry_name
                        )
                        self.create_sample_subgroup(
                            output_file,
                            raw_file_entry,
                            modified_entry_name,
                            sx_val,
                            sy_val,
                            sz_val,
                        )
                        self.create_entry_datasets(
                            output_file, raw_file_entry, modified_entry_name
                        )

    def main(self):
        """
        Main processing method for handling NX stress analysis on translation, radial,
        and potentially horizontal datasets based on the specified scan number for rotation.
        The method prepares output files, processes entries from the raw file,
        and handles different data processing modes based on the presence of rotational data.

        This method organizes data processing into translational and radial outputs by default,
        and includes horizontal outputs if a rotation scan number is specified.
        """
        directory, base_name = os.path.split(self.file_path)
        translational_output_file_name = os.path.join(
            directory, f"translational_NXstress_{base_name}"
        )
        radial_output_file_name = os.path.join(
            directory, f"radial_NXstress_{base_name}"
        )
        if self.scanNbForRotation is not None:
            horizontal_output_file_name = os.path.join(
                directory, f"horizontal_NXstress_{base_name}"
            )

        if os.path.exists(translational_output_file_name):
            os.remove(translational_output_file_name)
        if os.path.exists(radial_output_file_name):
            os.remove(radial_output_file_name)
        if self.scanNbForRotation is not None:
            if os.path.exists(horizontal_output_file_name):
                os.remove(horizontal_output_file_name)

        with h5py.File(self.file_path, "r") as raw_file:
            entries = list(raw_file)
            entries.sort(key=lambda x: float(x.rsplit(".", 2)[-2].rsplit("_", 1)[-1]))

            self.energy_calib = self.energy_calib_det0
            for entry_idx, entry in enumerate(entries):
                if (
                    self.scanNbForRotation is not None
                    and entry_idx + 1 >= self.scanNbForRotation
                ):
                    break
                self.process_scan(
                    translational_output_file_name, entry, detector_index=0
                )

            self.diffractogram_count = 0
            self.rangeFit = self.rangeFitDet1
            self.nbPeaksInBoxes = self.nbPeaksInBoxesDet1
            self.energy_calib = self.energy_calib_det1

            for entry_idx, entry in enumerate(entries):
                if (
                    self.scanNbForRotation is not None
                    and entry_idx + 1 == self.scanNbForRotation
                ):
                    break
                self.process_scan(radial_output_file_name, entry, detector_index=1)

            if self.scanNbForRotation is not None and self.rangeFitDet2 is not None:
                self.diffractogram_count = 0
                self.rangeFit = self.rangeFitDet2
                self.nbPeaksInBoxes = self.nbPeaksInBoxesDet2
                self.energy_calib = self.energy_calib_det0
                self.rotated = True
                for entry_idx, entry in enumerate(entries):
                    if entry_idx + 1 >= self.scanNbForRotation:
                        self.process_scan(
                            horizontal_output_file_name, entry, detector_index=0
                        )


class InteractivePlot:
    def __init__(self, data, title, start_index=0, total_length=None):
        self.data = data
        self.title = title
        self.fig, self.ax = plt.subplots()

        self.num_channels = len(data)
        self.start_index = start_index
        self.total_length = (
            total_length if total_length is not None else self.num_channels
        )
        self.energy_range = 300
        self.energy_per_channel = self.energy_range / self.total_length
        self.xs = [
            (channel + self.start_index) * self.energy_per_channel
            for channel in range(self.num_channels)
        ]
        self.ys = data

        self.ranges = []
        self.current_range = []

        # Immediately plot data
        self.plot_data()

    def plot_data(self):
        self.ax.plot(self.xs, self.ys, label="Data")
        self.ax.set_title(self.title)
        self.ax.set_xlabel("Energy (keV)")
        self.ax.set_ylabel("Counts")
        self.ax.legend()

    def onclick(self, event):
        if event.inaxes != self.ax:
            return
        if len(self.current_range) < 2:
            channel_index = int(
                (event.xdata - (self.start_index * self.energy_per_channel))
                / self.energy_per_channel
            )
            self.current_range.append(channel_index)
            self.ax.axvline(x=event.xdata, color="r", linestyle="--")
            self.fig.canvas.draw()

        if len(self.current_range) % 2 == 0:
            self.ranges.append(tuple(self.current_range))
            self.current_range = []

    def connect_events(self):
        self.cid = self.fig.canvas.mpl_connect("button_press_event", self.onclick)
        self.cid_escape = self.fig.canvas.mpl_connect(
            "key_press_event", self.on_key_press
        )

    def on_key_press(self, event):
        if event.key == "escape":
            self.fig.canvas.mpl_disconnect(self.cid)
            self.fig.canvas.mpl_disconnect(self.cid_escape)
            self.fig.canvas.draw_idle()

    def select_peak(self):
        # This method might be meant to prepare or finalize the interactive peak selection
        # Connect event handlers if not already connected
        if not hasattr(self, "cid"):
            self.cid = self.fig.canvas.mpl_connect("button_press_event", self.onclick)
            self.cid_escape = self.fig.canvas.mpl_connect(
                "key_press_event", self.on_key_press
            )

    def finalize_selection(self):
        # This can be called to disconnect event handlers and clean up
        if hasattr(self, "cid"):
            self.fig.canvas.mpl_disconnect(self.cid)
            self.fig.canvas.mpl_disconnect(self.cid_escape)
        return self.ranges


class GUI(tk.Tk):
    def __init__(self, data, handle_gui_data, title, prompt_hkl=False):
        super().__init__()
        self.data = data
        self.handle_gui_data = handle_gui_data
        self.geometry("1400x700")
        self.prompt_hkl = prompt_hkl
        self.title = title
        self.setup_ui()

    def setup_ui(self):
        standard_font = tkFont.Font(family="Helvetica", size=12)

        style = ttk.Style()
        style.configure("TLabel", font=standard_font)
        style.configure("TEntry", font=standard_font, padding=5)
        style.configure("TButton", font=standard_font, padding=5)

        if not self.prompt_hkl:
            self.plot = InteractivePlot(data=self.data, title=self.title)
            self.canvas = FigureCanvasTkAgg(self.plot.fig, master=self)
            self.canvas.draw()
            self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            self.plot.connect_events()
            toolbar = NavigationToolbar2Tk(self.canvas, self)
            toolbar.update()
            self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

            self.label = ttk.Label(
                self, text="Enter the number of peaks:", style="TLabel"
            )
            self.label.pack(padx=10, pady=10)
            self.peak_entry = ttk.Entry(self, style="TEntry")
            self.peak_entry.pack(padx=10, pady=10)
            self.peak_entry.bind(
                "<Return>", self.submit_peaks
            )  # Bind Enter key to submit_peaks

            self.submit_button = ttk.Button(
                self, text="Submit", style="TButton", command=self.submit_peaks
            )
            self.submit_button.pack(padx=10, pady=10)

            self.result_label = ttk.Label(self, text="", style="TLabel")
            self.result_label.pack(padx=10, pady=10)
        else:
            self.label_hkl = ttk.Label(
                self, text="Enter hkl value(s) separated by commas (e.g., '311, 211'):"
            )
            self.label_hkl.pack(padx=10, pady=10)
            self.hkl_entry = ttk.Entry(self)
            self.hkl_entry.pack(padx=10, pady=10)
            self.hkl_entry.bind("<Return>", self.submit_hkl)
            self.submit_hkl_button = ttk.Button(
                self, text="Submit hkl", command=self.submit_hkl
            )
            self.submit_hkl_button.pack(padx=10, pady=10)

    def submit_hkl(self, event=None):
        hkl_str = self.hkl_entry.get()
        hkl_values = hkl_str.split(",")
        valid = True
        h_miller, k_miller, l_miller = [], [], []
        for value in hkl_values:
            if len(value.strip()) == 3 and value.strip().isdigit():
                h_miller.append(int(value.strip()[0]))
                k_miller.append(int(value.strip()[1]))
                l_miller.append(int(value.strip()[2]))
            else:
                print("Invalid input:", value)
                valid = False
                break
        if valid:
            self.handle_gui_data(h_miller, k_miller, l_miller)
            self.quit()
            self.destroy()

    def submit_peaks(self, event=None):
        peaks = self.peak_entry.get()
        peak_ranges = self.plot.finalize_selection()
        self.result_label.config(text=f"Number of peaks: {peaks}")
        self.handle_gui_data(peak_ranges, peaks)
        self.quit()
        self.destroy()


def _read_config(file_path):
    with open(file_path, "r") as file:
        config = yaml.safe_load(file)
    return config


def initialize_system(config_path):
    config_values = _read_config(config_path)
    nx_stress = NXstressFromRaw(
        config_values["file_path"],
        det_calib_file_angle=config_values["det_calib_file_angle"],
        det_calib_file_energy=config_values["det_calib_file_energy"],
        with_cradle=config_values["with_cradle"],
        lattice=config_values["lattice"],
        phase_name=config_values["phase_name"],
        scanNbForRotation=int(config_values.get("scanNbForRotation", "None") or "0"),
        experimental_identifier=config_values["experimental_identifier"],
        collection_identifier=config_values["collection_identifier"],
    )
    return nx_stress
