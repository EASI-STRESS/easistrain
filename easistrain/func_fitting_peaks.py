import numpy as np
from func_fitting_functions import *
import scipy.optimize
import h5py


### This function fit a hkl peak saved in a result h5file and save the results in the same h5 file
### The results if the fit are saved in the 'fitting' group
### root_data: the path of the file where the h5 file is saved
### h5file: the name of the h5 file where the results are saved
### tth_min: 2theta start
### tth_max: 2theta end
### bgd_left: number of point to take at left (at the beginning of the tth range) for background fitting
### bgd_right: number of point to take at right (at the end of the tth range) for background fitting
### fct: the name of function to use for peak fitting: PVII, gauss, Lorentz or Ps
### hkl: crystallographic plane to fit
def fit(root_data, h5file, scan, tth_min, tth_max, bgd_left, bgd_right, fct, hkl, Rp):
    fh5_save = h5py.File(
        root_data + "/" + "Results" + "_" + h5file, "a"
    )  ### Create or append (if created) the file in which will be saved the results (integration, ...)
    r_h5file = h5py.File(
        root_data + "/" + "Results" + "_" + h5file, "r"
    )  ### Read the h5 file
    r_groups_scan = list(
        r_h5file.keys()
    )  # getting the list of the all the names of the scans in the h5 file
    if scan != "all":  # this loop executes if we define the name of the scan to process
        print("#*#*#*#*#*# Processing of the scan: " + scan + " #*#*#*#*#*#*")
        r_groups_images = list(
            r_h5file["/" + scan + "/raw_integration_2D"].keys()
        )  # getting the list images in a scan
        fit_group = fh5_save["/" + scan].create_group(
            "fitting_HKL=" + "(" + hkl + ")"
        )  # Creating of the fitting group (where the results of the fit will be saved)
        for nimg in range(
            1, np.shape(r_groups_images)[0]
        ):  # Iteration on the images in a scan
            print(
                "#*#*#*#*#*# Processing of the image: "
                + r_groups_images[nimg]
                + " #*#*#*#*#*#*"
            )
            print("######### Fitting started ############")
            img_fit = fit_group.create_group(
                r_groups_images[nimg]
            )  # creates a group 'image+nb' for each image
            tth_group = r_h5file[
                "/" + scan + "/raw_integration_2D/" + r_groups_images[nimg] + "/tth"
            ]  # read tth dataset matrix
            chi_group = r_h5file[
                "/" + scan + "/raw_integration_2D/" + r_groups_images[nimg] + "/chi"
            ]  # read chi dataset matrix
            cts_vs_tth_group = r_h5file[
                "/"
                + scan
                + "/raw_integration_2D/"
                + r_groups_images[nimg]
                + "/tth_vs_cts"
            ]  # read of the Intensity as function of two theta dataset matrix
            fit_matrix = np.zeros(
                (np.shape(chi_group)[0], 10), float
            )  # matrix on which will be stocked the results of the fit
            step = (tth_group[-1] - tth_group[0]) / (
                np.shape(tth_group)[0] - 1
            )  # The two theta step
            tth_range = tth_group[
                int(
                    np.round(((np.float(tth_min) / step) - (tth_group[0] / step)))
                ) : int(np.round((np.float(tth_max) / step) - (tth_group[0] / step)))
            ]  # definition of the two theta range as defined by tth_min and tth_max
            data_raw = cts_vs_tth_group[
                int(
                    np.round(((np.float(tth_min) / step) - (tth_group[0] / step)))
                ) : int(np.round((np.float(tth_max) / step) - (tth_group[0] / step))),
                :,
            ]  # the raw data for the range defined by tth_min and tth_max
            img_fit.create_dataset(
                "tth_range", dtype="f", data=tth_range
            )  # saving the tth range
            img_fit.create_dataset(
                "chi", dtype="f", data=chi_group
            )  # saving the chi (azim) angles
            data_bgd_subs = np.zeros(
                (np.shape(tth_range)[0], np.shape(chi_group)[0] + 1), float
            )  # the data for the range defined by tth_min and tth_max with bgd substracted
            data_fitted = np.zeros(
                (np.shape(tth_range)[0], np.shape(chi_group)[0] + 1), float
            )  # the fitted data
            data_error = np.zeros(
                (np.shape(tth_range)[0], np.shape(chi_group)[0] + 1), float
            )  # difference between fitted data and raw data
            tth_position = np.zeros(
                np.shape(chi_group)[0], float
            )  # matrix on which will be saved the tth positions
            for nazim in range(
                1, np.shape(cts_vs_tth_group)[1]
            ):  # iteration on the number of azimuth sectors (chi)
                # cts_range = data_raw[:,nazim] # Intensity corresponding to the two theta defined range
                # print(cts_range)
                # print(np.shape(data_raw))
                x_bgd = np.append(
                    data_raw[0 : int(bgd_left), 0],
                    data_raw[int(np.shape(data_raw)[0]) - int(bgd_right) : -1, 0],
                )  # tth points to be considered for bgd fitting
                y_bgd = np.append(
                    data_raw[0 : int(bgd_left), nazim],
                    data_raw[int(np.shape(data_raw)[0]) - int(bgd_right) : -1, nazim],
                )  # intensity points to be considered for bgd fitting
                # print(x_bgd)
                # print(y_bgd)
                # print(data_raw[0:int(bgd_left),nazim])
                # print(data_raw[int(bgd_right):-1, nazim])
                # print(np.shape(y_bgd))
                coef_bgd = np.polyfit(
                    x_bgd, y_bgd, 1
                )  # background fitting with a one degree polynome (a*x + b)
                bgd_func = np.poly1d(coef_bgd)  # the polynom used to fit the background
                fit_matrix[nazim - 1, 0] = coef_bgd[1]  # the constant b
                fit_matrix[nazim - 1, 1] = coef_bgd[0]  # the slope a
                data_bgd_subs[:, 0] = tth_range[:]  # tth
                data_bgd_subs[:, nazim] = data_raw[:, nazim] - bgd_func(
                    data_raw[:, 0]
                )  # substraction of the background
                # print('######### background substracted ############')
                #### definition of the first solution parameters (the vector P) of the fitting function
                P = np.zeros(4, float)
                P[0] = np.amax(data_bgd_subs[:, nazim])  # Max of intensity
                # print(P[0])
                # print(np.where(data_bgd_subs[:, nazim] == P[0]))
                P[1] = data_bgd_subs[
                    np.where(data_bgd_subs[:, nazim] == P[0])[0][0], 0
                ]  # tth
                # print(P[1])
                # print(np.where(np.round(data_bgd_subs[:, nazim]/(P[0]/2))==1))
                # print(np.where(np.round(data_bgd_subs[:, nazim]/(P[0]/2))==1)[0][0])
                # print('here')
                if (
                    np.size(
                        np.where(np.round(data_bgd_subs[:, nazim] / (P[0] / 2)) == 1)
                    )
                    == 0
                    or P[0] == 0
                ):
                    P[2] = 0.01
                    # print(np.size(np.where(np.round(data_bgd_subs[:, nazim]/(P[0]/2))==1)))
                else:
                    P[2] = np.abs(
                        data_bgd_subs[
                            np.where(
                                np.round(data_bgd_subs[:, nazim] / (P[0] / 2)) == 1
                            )[0][0],
                            0,
                        ]
                        - P[1]
                    )  # FWHM
                    # print(P[2])
                P[3] = 0.5
                ######## FITTING ###########
                x1 = data_bgd_subs[:, 0]
                y1 = data_bgd_subs[:, nazim]
                if fct == "PVII":
                    fitfunc = lambda P, x: func_PearsonVII(P, x)
                    errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                    R1, success = scipy.optimize.leastsq(
                        errfunc, P, args=(x1, y1), maxfev=10000
                    )
                    finalfit = fitfunc(R1, x1)
                if fct == "gauss":
                    fitfunc = lambda P, x: func_Gauss(P, x)
                    errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                    R1, success = scipy.optimize.leastsq(
                        errfunc, P, args=(x1, y1), maxfev=10000
                    )
                    finalfit = fitfunc(R1, x1)
                if fct == "Lorentz":
                    fitfunc = lambda P, x: func_Lorentz(P, x)
                    errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                    R1, success = scipy.optimize.leastsq(
                        errfunc, P, args=(x1, y1), maxfev=10000
                    )
                    finalfit = fitfunc(R1, x1)
                if fct == "PsV":
                    fitfunc = lambda P, x: func_pseudo_voigt(P, x)
                    errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                    R1, success = scipy.optimize.leastsq(
                        errfunc, P, args=(x1, y1), maxfev=10000
                    )
                    finalfit = fitfunc(R1, x1)
                else:
                    print("#### The fitting function is not defined")
                    print("#### please give fct as: PVII, gauss, Lorentz or PsV")
                    return
                fit_matrix[nazim - 1, 2] = R1[0]  # maximum peak intensity
                fit_matrix[nazim - 1, 3] = R1[1]  # tth position
                fit_matrix[nazim - 1, 4] = R1[2]  # FWHM or integral breadth
                fit_matrix[nazim - 1, 5] = R1[3]  # shape factor
                fit_matrix[nazim - 1, 6] = success  # fit success factor
                fit_matrix[nazim - 1, 8] = chi_group[nazim - 1]  # azimuth angle
                tth_position[nazim - 1] = R1[1]  # tth position
                data_fitted[:, 0] = tth_range[:]
                data_error[:, 0] = tth_range[:]
                data_fitted[:, nazim] = fitfunc(R1, x1[:]) + bgd_func(
                    data_raw[:, 0]
                )  # matrix containng the fitted data
                data_error[:, nazim] = np.abs(
                    data_raw[:, nazim] - data_fitted[:, nazim]
                )  # matrix containng absolute difference between fitted and raw data
                fit_matrix[nazim - 1, 7] = (
                    100 * np.sum(data_error[:, nazim]) / np.sum(data_fitted[:, nazim])
                )  # the Rp factor (goodness of the fit in %)
                fit_matrix[nazim - 1, 9] = np.sum(data_bgd_subs[:, nazim])
            img_fit.create_dataset("tth_position", dtype="f", data=tth_position)
            img_fit.create_dataset(
                "fit", dtype="f", data=fit_matrix
            )  # saving the chi (azim) angles
            img_fit.create_dataset(
                "data_bgd_subs", dtype="f", data=data_bgd_subs
            )  # saving data with bgd substracted
            img_fit.create_dataset(
                "data_fitted", dtype="f", data=data_fitted
            )  # saving fitted data
            img_fit.create_dataset(
                "data_raw", dtype="f", data=data_raw
            )  # saving raw data
            img_fit.create_dataset(
                "data_error", dtype="f", data=data_error
            )  # saving error between fitted and raw data
            print("######### Fitting completed ############")
            print("######### data saved in h5 file ############")
    else:  # This loop will be executed if the scan input is defined as all or whatever different from a name of scans in the h5 file (it processes all the scans in the h5 file)
        for iscan in range(np.shape(r_groups_scan)[0]):  # Iteration on the scans
            print(
                "#*#*#*#*#*# Processing of the scan: "
                + r_groups_scan[iscan]
                + " #*#*#*#*#*#*"
            )
            r_groups_images = list(
                r_h5file[r_groups_scan[iscan] + "/raw_integration_2D"].keys()
            )  # getting the list images in a scan
            fit_group = fh5_save[r_groups_scan[iscan]].create_group(
                "fitting_HKL=" + "(" + hkl + ")"
            )  # Creating of the fitting group (where the results of the fit will be saved)
            for nimg in range(
                1, np.shape(r_groups_images)[0]
            ):  # Iteration on the images in a scan
                print(
                    "#*#*#*#*#*# Processing of the image: "
                    + r_groups_images[nimg]
                    + " #*#*#*#*#*#*"
                )
                print("######### Fitting started ############")
                img_fit = fit_group.create_group(
                    r_groups_images[nimg]
                )  # creates a group 'image+nb' for each image
                tth_group = r_h5file[
                    r_groups_scan[iscan]
                    + "/raw_integration_2D/"
                    + r_groups_images[nimg]
                    + "/tth"
                ]  # read tth dataset matrix
                chi_group = r_h5file[
                    r_groups_scan[iscan]
                    + "/raw_integration_2D/"
                    + r_groups_images[nimg]
                    + "/chi"
                ]  # read chi dataset matrix
                cts_vs_tth_group = r_h5file[
                    r_groups_scan[iscan]
                    + "/raw_integration_2D/"
                    + r_groups_images[nimg]
                    + "/tth_vs_cts"
                ]  # read of the Intensity as function of two theta dataset matrix
                fit_matrix = np.zeros(
                    (np.shape(chi_group)[0], 10), float
                )  # matrix on which will be stocked the results of the fit
                step = (tth_group[-1] - tth_group[0]) / (
                    np.shape(tth_group)[0] - 1
                )  # The two theta step
                tth_range = tth_group[
                    int(
                        np.round(((np.float(tth_min) / step) - (tth_group[0] / step)))
                    ) : int(
                        np.round((np.float(tth_max) / step) - (tth_group[0] / step))
                    )
                ]  # definition of the two theta range as defined by tth_min and tth_max
                data_raw = cts_vs_tth_group[
                    int(
                        np.round(((np.float(tth_min) / step) - (tth_group[0] / step)))
                    ) : int(
                        np.round((np.float(tth_max) / step) - (tth_group[0] / step))
                    ),
                    :,
                ]  # the raw data for the range defined by tth_min and tth_max
                img_fit.create_dataset(
                    "tth_range", dtype="f", data=tth_range
                )  # saving the tth range
                img_fit.create_dataset(
                    "chi", dtype="f", data=chi_group
                )  # saving the chi (azim) angles
                data_bgd_subs = np.zeros(
                    (np.shape(tth_range)[0], np.shape(chi_group)[0] + 1), float
                )  # the data for the range defined by tth_min and tth_max with bgd substracted
                data_fitted = np.zeros(
                    (np.shape(tth_range)[0], np.shape(chi_group)[0] + 1), float
                )  # the fitted data
                data_error = np.zeros(
                    (np.shape(tth_range)[0], np.shape(chi_group)[0] + 1), float
                )  # difference between fitted data and raw data
                tth_position = np.zeros(
                    np.shape(chi_group)[0], float
                )  # matrix on which will be saved the tth positions
                for nazim in range(
                    1, np.shape(cts_vs_tth_group)[1]
                ):  # iteration on the number of azimuth sectors (chi)
                    # cts_range = data_raw[:,nazim] # Intensity corresponding to the two theta defined range
                    # print(cts_range)
                    # print(np.shape(data_raw))
                    x_bgd = np.append(
                        data_raw[0 : int(bgd_left), 0],
                        data_raw[int(np.shape(data_raw)[0]) - int(bgd_right) : -1, 0],
                    )  # tth points to be considered for bgd fitting
                    y_bgd = np.append(
                        data_raw[0 : int(bgd_left), nazim],
                        data_raw[
                            int(np.shape(data_raw)[0]) - int(bgd_right) : -1, nazim
                        ],
                    )  # intensity points to be considered for bgd fitting
                    # print(x_bgd)
                    # print(y_bgd)
                    # print(data_raw[0:int(bgd_left),nazim])
                    # print(data_raw[int(bgd_right):-1, nazim])
                    # print(np.shape(y_bgd))
                    coef_bgd = np.polyfit(
                        x_bgd, y_bgd, 1
                    )  # background fitting with a one degree polynome (a*x + b)
                    bgd_func = np.poly1d(
                        coef_bgd
                    )  # the polynom used to fit the background
                    fit_matrix[nazim - 1, 0] = coef_bgd[1]  # the constant b
                    fit_matrix[nazim - 1, 1] = coef_bgd[0]  # the slope a
                    data_bgd_subs[:, 0] = tth_range[:]  # tth
                    data_bgd_subs[:, nazim] = data_raw[:, nazim] - bgd_func(
                        data_raw[:, 0]
                    )  # substraction of the background
                    # print('######### background substracted ############')
                    #### definition of the first solution parameters (the vector P) of the fitting function
                    P = np.zeros(4, float)
                    P[0] = np.amax(data_bgd_subs[:, nazim])  # Max of intensity
                    # print(P[0])
                    # print(np.where(data_bgd_subs[:, nazim] == P[0]))
                    P[1] = data_bgd_subs[
                        np.where(data_bgd_subs[:, nazim] == P[0])[0][0], 0
                    ]  # tth
                    # print(P[1])
                    # print(np.where(np.round(data_bgd_subs[:, nazim]/(P[0]/2))==1))
                    # print(np.where(np.round(data_bgd_subs[:, nazim]/(P[0]/2))==1)[0][0])
                    # print('here')
                    if (
                        np.size(
                            np.where(
                                np.round(data_bgd_subs[:, nazim] / (P[0] / 2)) == 1
                            )
                        )
                        == 0
                        or P[0] == 0
                    ):
                        P[2] = 0.01
                        # print(np.size(np.where(np.round(data_bgd_subs[:, nazim]/(P[0]/2))==1)))
                    else:
                        P[2] = np.abs(
                            data_bgd_subs[
                                np.where(
                                    np.round(data_bgd_subs[:, nazim] / (P[0] / 2)) == 1
                                )[0][0],
                                0,
                            ]
                            - P[1]
                        )  # FWHM
                        # print(P[2])
                    P[3] = 0.5
                    ######## FITTING ###########
                    x1 = data_bgd_subs[:, 0]
                    y1 = data_bgd_subs[:, nazim]
                    if fct == "PVII":
                        fitfunc = lambda P, x: func_PearsonVII(P, x)
                        errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                        R1, success = scipy.optimize.leastsq(
                            errfunc, P, args=(x1, y1), maxfev=10000
                        )
                        finalfit = fitfunc(R1, x1)
                    if fct == "gauss":
                        fitfunc = lambda P, x: func_Gauss(P, x)
                        errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                        R1, success = scipy.optimize.leastsq(
                            errfunc, P, args=(x1, y1), maxfev=10000
                        )
                        finalfit = fitfunc(R1, x1)
                    if fct == "Lorentz":
                        fitfunc = lambda P, x: func_Lorentz(P, x)
                        errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                        R1, success = scipy.optimize.leastsq(
                            errfunc, P, args=(x1, y1), maxfev=10000
                        )
                        finalfit = fitfunc(R1, x1)
                    if fct == "PsV":
                        fitfunc = lambda P, x: func_pseudo_voigt(P, x)
                        errfunc = lambda P, x1, y1: fitfunc(P, x1) - y1
                        R1, success = scipy.optimize.leastsq(
                            errfunc, P, args=(x1, y1), maxfev=10000
                        )
                        finalfit = fitfunc(R1, x1)
                    else:
                        print("#### The fitting function is not defined")
                        print("#### please give fct as: PVII, gauss, Lorentz or PsV")
                        return
                    fit_matrix[nazim - 1, 2] = R1[0]  # maximum peak intensity
                    fit_matrix[nazim - 1, 3] = R1[1]  # tth position
                    fit_matrix[nazim - 1, 4] = R1[2]  # FWHM or integral breadth
                    fit_matrix[nazim - 1, 5] = R1[3]  # shape factor
                    fit_matrix[nazim - 1, 6] = success  # fit success factor
                    fit_matrix[nazim - 1, 8] = chi_group[nazim - 1]  # azimuth angle
                    tth_position[nazim - 1] = R1[1]  # tth position
                    data_fitted[:, 0] = tth_range[:]
                    data_error[:, 0] = tth_range[:]
                    data_fitted[:, nazim] = fitfunc(R1, x1[:]) + bgd_func(
                        data_raw[:, 0]
                    )  # matrix containng the fitted data
                    data_error[:, nazim] = np.abs(
                        data_raw[:, nazim] - data_fitted[:, nazim]
                    )  # matrix containng absolute difference between fitted and raw data
                    fit_matrix[nazim - 1, 7] = (
                        100
                        * np.sum(data_error[:, nazim])
                        / np.sum(data_fitted[:, nazim])
                    )  # the Rp factor (goodness of the fit in %)
                    fit_matrix[nazim - 1, 9] = np.sum(data_bgd_subs[:, nazim])
                img_fit.create_dataset("tth_position", dtype="f", data=tth_position)
                img_fit.create_dataset(
                    "fit", dtype="f", data=fit_matrix
                )  # saving the chi (azim) angles
                img_fit.create_dataset(
                    "data_bgd_subs", dtype="f", data=data_bgd_subs
                )  # saving data with bgd substracted
                img_fit.create_dataset(
                    "data_fitted", dtype="f", data=data_fitted
                )  # saving fitted data
                img_fit.create_dataset(
                    "data_raw", dtype="f", data=data_raw
                )  # saving raw data
                img_fit.create_dataset(
                    "data_error", dtype="f", data=data_error
                )  # saving error between fitted and raw data
                print("######### Fitting completed ############")
                print("######### data saved in h5 file ############")
    return


def clean_fit(
    root_data, h5file, scan, tth_min, tth_max, bgd_left, bgd_right, fct, hkl, Rp
):
    fh5_save = h5py.File(
        root_data + "/" + "Results" + "_" + h5file, "a"
    )  ### Read the h5 file
    r_h5file = h5py.File(
        root_data + "/" + "Results" + "_" + h5file, "r"
    )  ### Create the file in which will be saved the results (integration, ...)
    r_groups_scan = list(
        r_h5file.keys()
    )  # getting the list of the all the names of the scans in the h5 file
    if scan != "all":
        print("#*#*#*#*#*# Processing of the scan: " + scan + " #*#*#*#*#*#*")
        r_groups_images = list(
            r_h5file["/" + scan + "/fitting_HKL=" + "(" + hkl + ")"].keys()
        )  # getting the list images in a scan
        clean_group = fh5_save["/" + scan].create_group(
            "fitting_HKL=" + "(" + hkl + ")_cleaned"
        )  # Creating of the group where the cleaned results will be saved
        for nimg in range(
            np.shape(r_groups_images)[0]
        ):  # Iteration on the images in a scan
            print(
                "#*#*#*#*#*# Processing of the image: "
                + r_groups_images[nimg]
                + " #*#*#*#*#*#*"
            )
            print("######### Cleaning started ############")
            img_fit_cleaned = clean_group.create_group(
                r_groups_images[nimg]
            )  # creates a group 'image+nb' for each image
            uncleaned_fit_matrix = r_h5file[
                "/"
                + scan
                + "/fitting_HKL="
                + "("
                + hkl
                + ")/"
                + r_groups_images[nimg]
                + "/fit"
            ]  # The fit matrix (where the results of the non-cleaned fit are saved)
            uncleaned_tth_position = r_h5file[
                "/"
                + scan
                + "/fitting_HKL="
                + "("
                + hkl
                + ")/"
                + r_groups_images[nimg]
                + "/tth_position"
            ]  # The fit matrix (where the results of the non-cleaned fit are saved)
            uncleaned_chi = r_h5file[
                "/"
                + scan
                + "/fitting_HKL="
                + "("
                + hkl
                + ")/"
                + r_groups_images[nimg]
                + "/chi"
            ]  # The fit matrix (where the results of the non-cleaned fit are saved)
            # cleaned_fit_matrix = np.zeros((np.shape(uncleaned_fit_matrix)[0],np.shape(uncleaned_fit_matrix)[1]),float) # Matrix on which will be stocked the cleaned fit matrix
            index = np.array([])
            for nazim in range(np.shape(uncleaned_fit_matrix)[0]):
                if fct == "PsV":
                    if (
                        any(uncleaned_fit_matrix[nazim, 2:8] <= 0)
                        or uncleaned_fit_matrix[nazim, 5] < 0
                        or uncleaned_fit_matrix[nazim, 5] > 1
                        or uncleaned_fit_matrix[nazim, 7] > int(Rp)
                        or uncleaned_fit_matrix[nazim, 2]
                        + (
                            uncleaned_fit_matrix[nazim, 0]
                            + (
                                uncleaned_fit_matrix[nazim, 1]
                                * uncleaned_fit_matrix[nazim, 3]
                            )
                        )
                        <= (
                            uncleaned_fit_matrix[nazim, 0]
                            + (
                                uncleaned_fit_matrix[nazim, 1]
                                * uncleaned_fit_matrix[nazim, 3]
                            )
                        )
                        + (
                            3
                            * np.sqrt(
                                (
                                    uncleaned_fit_matrix[nazim, 0]
                                    + (
                                        uncleaned_fit_matrix[nazim, 1]
                                        * uncleaned_fit_matrix[nazim, 3]
                                    )
                                )
                            )
                        )
                    ):
                        index = np.append(index, nazim)
                else:
                    if (
                        any(uncleaned_fit_matrix[nazim, 2:8] <= 0)
                        or uncleaned_fit_matrix[nazim, 7] > int(Rp)
                        or uncleaned_fit_matrix[nazim, 2]
                        + (
                            uncleaned_fit_matrix[nazim, 0]
                            + (
                                uncleaned_fit_matrix[nazim, 1]
                                * uncleaned_fit_matrix[nazim, 3]
                            )
                        )
                        <= (
                            uncleaned_fit_matrix[nazim, 0]
                            + (
                                uncleaned_fit_matrix[nazim, 1]
                                * uncleaned_fit_matrix[nazim, 3]
                            )
                        )
                        + (
                            3
                            * np.sqrt(
                                (
                                    uncleaned_fit_matrix[nazim, 0]
                                    + (
                                        uncleaned_fit_matrix[nazim, 1]
                                        * uncleaned_fit_matrix[nazim, 3]
                                    )
                                )
                            )
                        )
                    ):
                        index = np.append(index, nazim)
            img_fit_cleaned.create_dataset(
                "fit_cleaned",
                dtype="f",
                data=np.delete(uncleaned_fit_matrix, np.int32(index), 0),
            )
            img_fit_cleaned.create_dataset(
                "tth_position_cleaned",
                dtype="f",
                data=np.delete(uncleaned_tth_position, np.int32(index), 0),
            )
            img_fit_cleaned.create_dataset(
                "chi_cleaned",
                dtype="f",
                data=np.delete(uncleaned_chi, np.int32(index), 0),
            )
            print("######### Cleaning completed ############")
            print("######### data saved in h5 file ############")
    else:
        for iscan in range(np.shape(r_groups_scan)[0]):  # Iteration on the scans
            print(
                "#*#*#*#*#*# Processing of the scan: "
                + r_groups_scan[iscan]
                + " #*#*#*#*#*#*"
            )
            r_groups_images = list(
                r_h5file[
                    r_groups_scan[iscan] + "/fitting_HKL=" + "(" + hkl + ")"
                ].keys()
            )  # getting the list images in a scan
            clean_group = fh5_save[r_groups_scan[iscan]].create_group(
                "fitting_HKL=" + "(" + hkl + ")_cleaned"
            )  # Creating of the group where the cleaned results will be saved
            for nimg in range(
                np.shape(r_groups_images)[0]
            ):  # Iteration on the images in a scan
                print(
                    "#*#*#*#*#*# Processing of the image: "
                    + r_groups_images[nimg]
                    + " #*#*#*#*#*#*"
                )
                print("######### Cleaning started ############")
                img_fit_cleaned = clean_group.create_group(
                    r_groups_images[nimg]
                )  # creates a group 'image+nb' for each image
                uncleaned_fit_matrix = r_h5file[
                    r_groups_scan[iscan]
                    + "/fitting_HKL="
                    + "("
                    + hkl
                    + ")/"
                    + r_groups_images[nimg]
                    + "/fit"
                ]  # The fit matrix (where the results of the non-cleaned fit are saved)
                uncleaned_tth_position = r_h5file[
                    r_groups_scan[iscan]
                    + "/fitting_HKL="
                    + "("
                    + hkl
                    + ")/"
                    + r_groups_images[nimg]
                    + "/tth_position"
                ]  # The fit matrix (where the results of the non-cleaned fit are saved)
                uncleaned_chi = r_h5file[
                    r_groups_scan[iscan]
                    + "/fitting_HKL="
                    + "("
                    + hkl
                    + ")/"
                    + r_groups_images[nimg]
                    + "/chi"
                ]  # The fit matrix (where the results of the non-cleaned fit are saved)
                # cleaned_fit_matrix = np.zeros((np.shape(uncleaned_fit_matrix)[0],np.shape(uncleaned_fit_matrix)[1]),float) # Matrix on which will be stocked the cleaned fit matrix
                index = np.array([])
                for nazim in range(np.shape(uncleaned_fit_matrix)[0]):
                    if fct == "PsV":
                        if (
                            any(uncleaned_fit_matrix[nazim, 2:8] <= 0)
                            or uncleaned_fit_matrix[nazim, 5] < 0
                            or uncleaned_fit_matrix[nazim, 5] > 1
                            or uncleaned_fit_matrix[nazim, 7] > int(Rp)
                            or uncleaned_fit_matrix[nazim, 2]
                            + (
                                uncleaned_fit_matrix[nazim, 0]
                                + (
                                    uncleaned_fit_matrix[nazim, 1]
                                    * uncleaned_fit_matrix[nazim, 3]
                                )
                            )
                            <= (
                                uncleaned_fit_matrix[nazim, 0]
                                + (
                                    uncleaned_fit_matrix[nazim, 1]
                                    * uncleaned_fit_matrix[nazim, 3]
                                )
                            )
                            + (
                                3
                                * np.sqrt(
                                    (
                                        uncleaned_fit_matrix[nazim, 0]
                                        + (
                                            uncleaned_fit_matrix[nazim, 1]
                                            * uncleaned_fit_matrix[nazim, 3]
                                        )
                                    )
                                )
                            )
                        ):
                            index = np.append(index, nazim)
                    else:
                        if (
                            any(uncleaned_fit_matrix[nazim, 2:8] <= 0)
                            or uncleaned_fit_matrix[nazim, 7] > int(Rp)
                            or uncleaned_fit_matrix[nazim, 2]
                            + (
                                uncleaned_fit_matrix[nazim, 0]
                                + (
                                    uncleaned_fit_matrix[nazim, 1]
                                    * uncleaned_fit_matrix[nazim, 3]
                                )
                            )
                            <= (
                                uncleaned_fit_matrix[nazim, 0]
                                + (
                                    uncleaned_fit_matrix[nazim, 1]
                                    * uncleaned_fit_matrix[nazim, 3]
                                )
                            )
                            + (
                                3
                                * np.sqrt(
                                    (
                                        uncleaned_fit_matrix[nazim, 0]
                                        + (
                                            uncleaned_fit_matrix[nazim, 1]
                                            * uncleaned_fit_matrix[nazim, 3]
                                        )
                                    )
                                )
                            )
                        ):
                            index = np.append(index, nazim)
                img_fit_cleaned.create_dataset(
                    "fit_cleaned",
                    dtype="f",
                    data=np.delete(uncleaned_fit_matrix, np.int32(index), 0),
                )
                img_fit_cleaned.create_dataset(
                    "tth_position_cleaned",
                    dtype="f",
                    data=np.delete(uncleaned_tth_position, np.int32(index), 0),
                )
                img_fit_cleaned.create_dataset(
                    "chi_cleaned",
                    dtype="f",
                    data=np.delete(uncleaned_chi, np.int32(index), 0),
                )
                print("######### Cleaning completed ############")
                print("######### data saved in h5 file ############")
    return
