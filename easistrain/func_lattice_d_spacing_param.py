import pyFAI
import numpy as np
import h5py


def lattice_param(root_data, h5file, poni_file, h, k, l):
    ai = pyFAI.load(poni_file)  # loading the poni file (integration geometry)
    wlgth = (ai.wavelength) * 10 ** 10  # the wavelength in angstrom
    fh5_save = h5py.File(
        root_data + "/" + "Results" + "_" + h5file, "a"
    )  ### Create the file in which will be saved the results (integration, ...)
    r_h5file = h5py.File(
        root_data + "/" + "Results" + "_" + h5file, "r"
    )  ### Read the h5 file
    r_groups_scan = list(
        r_h5file.keys()
    )  # getting the list of all the names of the scans in the h5 file
    for iscan in range(np.shape(r_groups_scan)[0]):  # Iteration on the scans
        print(
            "#*#*#*#*#*# Processing of the scan: "
            + r_groups_scan[iscan]
            + " #*#*#*#*#*#*"
        )
        r_groups_images = list(
            r_h5file[
                r_groups_scan[iscan]
                + "/fitting_HKL="
                + "("
                + str(h)
                + str(k)
                + str(l)
                + ")_cleaned"
            ].keys()
        )  # getting the list images in a scan
        save1 = fh5_save[r_groups_scan[iscan]].create_group(
            "latt_param_d_spacing_" + "(" + str(h) + str(k) + str(l) + ")"
        )  # Creating of the group where the results will be saved
        for nimg in range(
            np.shape(r_groups_images)[0]
        ):  # Iteration on the images in a scan
            print(
                "#*#*#*#*#*# Processing of the image: "
                + r_groups_images[nimg]
                + " #*#*#*#*#*#*"
            )
            print(
                "######### Calculation of lattice parameter and d-spacing started ############"
            )
            img_save1 = save1.create_group(
                r_groups_images[nimg]
            )  # creates a group 'image+nb' for each image
            tth_pos = r_h5file[
                r_groups_scan[iscan]
                + "/fitting_HKL="
                + "("
                + str(h)
                + str(k)
                + str(l)
                + ")_cleaned/"
                + r_groups_images[nimg]
                + "/tth_position_cleaned"
            ]  # read tth_positions dataset matrix
            chi_cleaned = r_h5file[
                r_groups_scan[iscan]
                + "/fitting_HKL="
                + "("
                + str(h)
                + str(k)
                + str(l)
                + ")_cleaned/"
                + r_groups_images[nimg]
                + "/chi_cleaned"
            ]  # read tth_positions dataset matrix
            a_matrix = np.zeros(
                (np.shape(tth_pos)[0], 2), float
            )  # matrix to stock the calculated lattice parameters
            d_spacing_matrix = np.zeros(
                (np.shape(tth_pos)[0], 2), float
            )  # matrix to stock the calculated d_spacings
            d_spacing_matrix[:, 0] = wlgth / (
                2 * np.sin(np.radians(0.5 * tth_pos[:]))
            )  # Calculation of the d_spacing
            d_spacing_matrix[:, 1] = chi_cleaned[:]
            a_matrix[:, 0] = d_spacing_matrix[:, 0] * np.sqrt(
                (int(h) ** 2) + (int(k) ** 2) + (int(l) ** 2)
            )  # Calculation of the lattice parameter
            a_matrix[:, 1] = chi_cleaned[:]
            a_mean = np.mean(
                a_matrix[:, 0]
            )  # Calculation o the mean lattice parameter value
            a_std = np.std(
                a_matrix[:, 0]
            )  # Calculation of the std of the lattice parameter
            d_mean = np.mean(
                d_spacing_matrix[:, 0]
            )  # Calculation of the mean d-spacing value
            d_std = np.std(
                d_spacing_matrix[:, 0]
            )  # Calculation of the std of the d-spacing value
            print(
                "######### Calculation of lattice parameter and d-spacing completed ############"
            )
            img_save1.create_dataset("d_spacing", dtype="f", data=d_spacing_matrix)
            img_save1.create_dataset("d_spacing_mean", dtype="f", data=d_mean)
            img_save1.create_dataset("d_spacing_std", dtype="f", data=d_std)
            img_save1.create_dataset("latt_param", dtype="f", data=a_matrix)
            img_save1.create_dataset("latt_param_mean", dtype="f", data=a_mean)
            img_save1.create_dataset("latt_param_std", dtype="f", data=a_std)
            print("######### data saved in h5 file ############")
    return
