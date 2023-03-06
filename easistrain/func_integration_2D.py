import os
import numpy
import h5py
import fabio
import pyFAI
import hdf5plugin  # noqa


def integration_2D(
    root_data,
    h5file,
    scan,
    numScan,
    detector_name,
    poni_file,
    npt_rad,
    npt_azim,
    x_unit,
    im_dark,
    im_mask,
    rad_range,
    azim_range,
    errorModel,
    imFlat,
    gon1,
    gon2,
    gon3,
    chiGon1,
    omegaGon2,
    phiGon3,
):
    """Integrates a 2D image using integrate2D pyfai's function and save the results in a h5file named Results_name of the h5 file containing the image

    1) It looks on the image on the h5 file
    2) creates a mask based on the int_max and int_min the user gives
    3) Integrate the 2D image in number of given sectors (npt_azim) and with a given number of points in the radial direction (npt_rad)
    4) saves the results in a h5file named Results_name of the h5 file containing the image
    Definition of the function inputs:
    :param root_data: the path of the file where the h5 file is saved
    :param h5file: the name of the h5 file where the images are saved
    :param scan: The name of the group (in h5file) on which the concerned measurement are saved
    :param detector name: name of the used detector
    :param poni_file: the path of the poni file which define the experimental geometry (ai)
    :param int_min: minimum intensity below which the pixels have to be masked (deleted)
    :param int_max: maximum intensity above which the pixels have to be masked (deleted)
    :param npt_rad: number of points in the radial direction (number of 2theta points)
    :param npt_azim: number of sectors for the integration in the azimutha direction (around the ring)
    :param x_unit: the unity of the x axis (2th_deg or 2th_rad or q_A^-1...)
    :param im_dark: the name of the dark image (if no dark image exist please give 0 as argument)
    :param im_mask: the name of the mask image (if no mask image exist please give 0 as argument)
    :param rad_range: the radial range to integrate in the radial direction (2tth, q, d, ...)
    :param azim_range: the range of azimuthal range (if we want to integrate just a portion of DS ring)
    """
    print(im_dark)
    with h5py.File(
        os.path.join(root_data, "Results" + "_" + h5file), "a"
    ) as fh5_save:  ### Create th file in which will be saved the results (integration, ...)
        with h5py.File(
            os.path.join(root_data, h5file), "r"
        ) as r_h5file:  ## Read the h5 file of raw data
            r_groups_scan = list(
                r_h5file.keys()
            )  # getting the list of the all the names of the scans in the h5 file
            ai = pyFAI.load(
                poni_file
            )  ### Load the poni file describing the integration geometry
            with open(os.path.join(root_data, "exe_integration.log"), "a") as flog:
                # if scan != 'all': # this loop executes if we define the name of the scan to process
                # 	print('#*#*#*#*#*#*#*#*# Processing of the scan: ' + scan + ' #*#*#*#*#*#*#*#*#')
                # 	flog.write('#*#*#*#*#*# Processing of the scan: ' + scan + ' #*#*#*#*#*#*#*#*# \n')
                # 	level_1 = fh5_save.create_group(scan) ### Create the group on which will be putted the data of integration
                # 	level_1_subg_1 = level_1.create_group('raw_integration_2D')
                # 	level_1_subg_3 = level_1_subg_1.create_group('Integration_parameter')
                # 	level_1_subg_4 = level_1_subg_1.create_group('Position')
                # 	level_1_subg_3.create_dataset('npt_azim', dtype ='f', data = int(npt_azim))
                # 	level_1_subg_3.create_dataset('npt_rad', dtype ='f', data = int(npt_rad))
                # 	image_nc = r_h5file['/' + scan + '/measurement/' + detector_name] ### Read the image or images
                # 	image = numpy.float64(image_nc) # convert the image matrix from int32 to float
                # 	g1_axis = r_h5file['/' + scan + '/instrument/positioners/' + gon1][()] # position of the sample along the beam
                # 	g2_axis = r_h5file['/' + scan + '/instrument/positioners/' + gon2][()] # position of the samlple in the dir G2
                # 	g3_axis = r_h5file['/' + scan + '/instrument/positioners/' + gon3][()] # position of the samlple in the dir G3 (perpendicular to the surface of the sample)
                # 	rg1_chi = r_h5file['/' + scan + '/instrument/positioners/' + chiGon1][()] # rotation around G1
                # 	rg2_omega = r_h5file['/' + scan + '/instrument/positioners/' + omegaGon2][()] # rotation around G2
                # 	rg3_phi = r_h5file['/' + scan + '/instrument/positioners/' + phiGon3][()] # rotation around G2
                # 	print('### The matrix image is a ' +  str(numpy.ndim(image)) + 'D matrix')
                # 	flog.write('### The matrix image is a ' +  str(numpy.ndim(image)) + 'D matrix \n')
                # 	if (numpy.ndim(image) == 2): # do a test on the dimension of the image matrix
                # 		rslt_matrix_cts = numpy.zeros((int(npt_rad), int(npt_azim) + 1), float)
                # 		rslt_matrix_chi = numpy.zeros((int(npt_azim)), float)
                # 		rslt_matrix_tth = numpy.zeros((int(npt_rad)), float)
                # 		if im_mask is None:
                # 		    #mask_mat = numpy.zeros((numpy.shape(image)[0], numpy.shape(image)[1]), float) # creating a zeroes matrix as no mask will be used in this case
                # 		    print('### No mask was used for the integration')
                # 		    flog.write('### No mask was used for the integration \n')
                # 		else:
                # 		    mask_mat = fabio.open(im_mask).data
                # 		    print('### The image: ' + im_mask + ' ' + 'was used as mask')
                # 		    flog.write('### The image: ' + im_mask + ' ' + 'was used as mask \n')
                # 		#print('### Do not worry I am creating the mask')
                # 		# for ix in range(numpy.shape(image)[0]): # loop on the image dimension and generation of the mask
                # 		# 	for iy in range(numpy.shape(image)[1]):
                # 		# 		if (image[ix,iy] < numpy.float(int_min) or image[ix,iy] > numpy.float(int_max)):
                # 		# 			mask_mat[ix,iy] = 1
                # 		# 			#print(str(ix) + '/' + str(iy))
                # 		#print('### Hola, I finished creating the mask')
                # 		print('### Integration started with a'+ ' ' + npt_rad + ' ' + 'steps in tth and with a' + ' ' + npt_azim + ' ' + 'sectors')
                # 		flog.write('### Integration started with a'+ ' ' + npt_rad + ' ' + 'steps in tth and with a' + ' ' + npt_azim + ' ' + 'sectors \n')
                # 		cts, tth, chi = ai.integrate2d(image, int(npt_rad), int(npt_azim), correctSolidAngle = True, polarization_factor=0.95, unit=x_unit, method="splitpixel", mask=mask_mat, dark = im_dark, radial_range = rad_range, azimuth_range = azim_range, error_model = errorModel, flat = imFlat)
                # 		print('### Integration copmleted')
                # 		flog.write('### Integration copmleted \n')
                # 		#print(numpy.shape(tth))
                # 		#print(numpy.shape(cts))
                # 		#print(numpy.shape(rslt_matrix_cts))
                # 		rslt_matrix_cts[:,0] = tth[:] # putting the x axis variables (2tth or q ...) in the first column of the result matrix
                # 		rslt_matrix_chi[:] = chi[:] # putting the chi angles (around the DS ring) in the chi dataset matrix
                # 		rslt_matrix_tth[:] = tth[:] # putting the x axis variables (2tth or q ...) in its dataset matrix
                # 		print('### Saving in h5file started')
                # 		flog.write('### Saving in h5file started \n')
                # 		for ichi in range(numpy.shape(chi)[0]):
                # 			rslt_matrix_cts[:,ichi+1] = cts[ichi,:] # saving the binned intensity as function of chi angles
                # 		level_1_subg_2 = level_1_subg_1.create_group('image_00000')
                # 		level_1_subg_2.create_dataset('tth_vs_cts', dtype ='f', data = rslt_matrix_cts)
                # 		level_1_subg_2.create_dataset('chi', dtype ='f', data = rslt_matrix_chi)
                # 		level_1_subg_2.create_dataset('tth', dtype ='f', data = rslt_matrix_tth)
                # 		level_1_subg_4.create_dataset('S1', dtype ='f', data = -g1_axis)
                # 		level_1_subg_4.create_dataset('S2', dtype ='f', data = -g2_axis)
                # 		level_1_subg_4.create_dataset('S3', dtype ='f', data = -g3_axis)
                # 		level_1_subg_4.create_dataset('Chi', dtype ='f', data = rg1_chi)
                # 		level_1_subg_4.create_dataset('Omega', dtype ='f', data = rg2_omega)
                # 		level_1_subg_4.create_dataset('phi', dtype ='f', data = rg3_phi)
                # 		print('### Saving in h5file completed')
                # 		flog.write('### Saving in h5file completed \n')
                # 	else: # is executed if the matrix dimension is not 2 (meaninig if it contains more than one image)
                # 		#print('### The matrix image is a ' +  str(numpy.ndim(image)) + 'D matrix')
                # 		if im_mask is None:
                # 		    #mask_mat = numpy.zeros((numpy.shape(image)[0], numpy.shape(image)[1]), float) # creating a zeroes matrix to be filled later with the mask
                # 		    print('### No mask was used for the integration')
                # 		    flog.write('### No mask was used for the integration \n')
                # 		else:
                # 		    mask_mat = fabio.open(im_mask).data
                # 		    print('### The image: ' + im_mask + ' ' + 'was used as mask')
                # 		    flog.write('### The image: ' + im_mask + ' ' + 'was used as mask \n')
                # 		for i in range(numpy.shape(image)[0]):
                # 			print('*#*#*#*#* Processing of image_' + '%.0f'%i + '*#*#*#*#*')
                # 			flog.write('*#*#*#*#* Processing of image_' + '%.0f'%i + '*#*#*#*#* \n')
                # 			rslt_matrix_cts = numpy.zeros((int(npt_rad), int(npt_azim) + 1), float)
                # 			rslt_matrix_chi = numpy.zeros((int(npt_azim)), float)
                # 			rslt_matrix_tth = numpy.zeros((int(npt_rad)), float)
                # 			#mask_mat = numpy.zeros((numpy.shape(image)[0], numpy.shape(image)[1]), float)
                # 			#print('### Do not worry I am creating the mask')
                # 			# for ix in range(numpy.shape(image)[1]):
                # 			# 	for iy in range(numpy.shape(image)[2]):
                # 			# 		if (image[i,ix,iy] < numpy.float(int_min) or image[i,ix,iy] > numpy.float(int_max)):
                # 			# 			mask_mat[ix,iy] = 1
                # 			#print('### Hola, I finished creating the mask')
                # 			print('### Integration started with a' + ' ' + npt_rad + ' ' + 'steps in tth and with a' + ' ' + npt_azim + ' ' + 'sectors')
                # 			flog.write('### Integration started with a' + ' ' + npt_rad + ' ' + 'steps in tth and with a' + ' ' + npt_azim + ' ' + 'sectors \n')
                # 			cts, tth, chi = ai.integrate2d(image[i], int(npt_rad), int(npt_azim), correctSolidAngle = True, polarization_factor=0.95, unit = x_unit, method = 'splitpixel', mask = mask_mat, dark = im_dark, radial_range = rad_range, azimuth_range = azim_range, error_model = errorModel, flat = imFlat)
                # 			print('### Integration copmleted')
                # 			flog.write('### Integration copmleted \n')
                # 			rslt_matrix_cts[:,0] = tth[:]
                # 			rslt_matrix_chi[:] = chi[:]
                # 			rslt_matrix_tth[:] = tth[:]
                # 			print('### Saving in h5file started')
                # 			flog.write('### Saving in h5file started \n')
                # 			for ichi in range(numpy.shape(chi)[0]):
                # 				rslt_matrix_cts[:,ichi+1] = cts[ichi,:]
                # 			if (i < 10):
                # 				level_1_subg_2 = level_1_subg_1.create_group('image' + '_0000' + '%.0f'%i)
                # 			if ((i>9) and (i<100)):
                # 				level_1_subg_2 = level_1_subg_1.create_group('image' + '_000' + '%.0f'%i)
                # 			if ((i>99) and (i<1000)):
                # 				level_1_subg_2 = level_1_subg_1.create_group('image' + '_00' + '%.0f'%i)
                # 			if ((i>999) and (i<10000)):
                # 				level_1_subg_2 = level_1_subg_1.create_group('image' + '_0' + '%.0f'%i)
                # 			if ((i>9999) and (i<100000)):
                # 				level_1_subg_2 = level_1_subg_1.create_group('image' + '_' + '%.0f'%i)
                # 			level_1_subg_2.create_dataset('tth_vs_cts', dtype ='f', data = rslt_matrix_cts)
                # 			level_1_subg_2.create_dataset('chi', dtype ='f', data = rslt_matrix_chi)
                # 			level_1_subg_2.create_dataset('tth', dtype ='f', data = rslt_matrix_tth)
                # 			print('### Saving in h5file completed')
                # 			flog.write('### Saving in h5file completed \n')
                # 		level_1_subg_4.create_dataset('S1', dtype ='f', data = -g1_axis)
                # 		level_1_subg_4.create_dataset('S2', dtype ='f', data = -g2_axis)
                # 		level_1_subg_4.create_dataset('S3', dtype ='f', data = -g3_axis)
                # 		level_1_subg_4.create_dataset('Chi', dtype ='f', data = rg1_chi)
                # 		level_1_subg_4.create_dataset('Omega', dtype ='f', data = rg2_omega)
                # 		level_1_subg_4.create_dataset('phi', dtype ='f', data = rg3_phi)
                if (
                    scan == "all" and numScan is None
                ):  # This loop will be executed if the scan=all and numScan=None (it processes all the scans in the h5 file)
                    for iscan in range(
                        numpy.shape(r_groups_scan)[0]
                    ):  # Iteration on the scans
                        print(
                            "#*#*#*#*#*# Processing of the scan: "
                            + r_groups_scan[iscan]
                            + " #*#*#*#*#*#*"
                        )
                        flog.write(
                            "#*#*#*#*#*# Processing of the scan: "
                            + r_groups_scan[iscan]
                            + " #*#*#*#*#*#\n"
                        )
                        if (
                            "measurement" in r_h5file[r_groups_scan[iscan]].keys()
                            and detector_name
                            in r_h5file[r_groups_scan[iscan] + "/measurement"].keys()
                        ):
                            level_1 = fh5_save.create_group(
                                r_groups_scan[iscan]
                            )  ### Create the group on which will be putted the data of integration
                            level_1_subg_1 = level_1.create_group("raw_integration_2D")
                            level_1_subg_3 = level_1_subg_1.create_group(
                                "Integration_parameter"
                            )
                            level_1_subg_3.create_dataset(
                                "npt_azim", dtype="f", data=int(npt_azim)
                            )
                            level_1_subg_3.create_dataset(
                                "npt_rad", dtype="f", data=int(npt_rad)
                            )
                            image_nc = r_h5file[
                                "/"
                                + r_groups_scan[iscan]
                                + "/measurement/"
                                + detector_name
                            ]  ### Read the image or images
                            image = numpy.float64(
                                image_nc
                            )  # convert the image matrix from int32 to float
                            if gon1 is not None:
                                g1_axis = r_h5file[
                                    "/"
                                    + r_groups_scan[iscan]
                                    + "/instrument/positioners/"
                                    + gon1
                                ][
                                    ()
                                ]  # position of the sample along the beam
                            if gon2 is not None:
                                g2_axis = r_h5file[
                                    "/"
                                    + r_groups_scan[iscan]
                                    + "/instrument/positioners/"
                                    + gon2
                                ][
                                    ()
                                ]  # position of the samlple in the dir G2
                            if gon3 is not None:
                                g3_axis = r_h5file[
                                    "/"
                                    + r_groups_scan[iscan]
                                    + "/instrument/positioners/"
                                    + gon3
                                ][
                                    ()
                                ]  # position of the samlple in the dir G3 (perpendicular to the surface of the sample)
                            if chiGon1 is not None:
                                rg1_chi = r_h5file[
                                    "/"
                                    + r_groups_scan[iscan]
                                    + "/instrument/positioners/"
                                    + chiGon1
                                ][
                                    ()
                                ]  # rotation around G1
                            if omegaGon2 is not None:
                                rg2_omega = r_h5file[
                                    "/"
                                    + r_groups_scan[iscan]
                                    + "/instrument/positioners/"
                                    + omegaGon2
                                ][
                                    ()
                                ]  # rotation around G2
                            if phiGon3 is not None:
                                rg3_phi = r_h5file[
                                    "/"
                                    + r_groups_scan[iscan]
                                    + "/instrument/positioners/"
                                    + phiGon3
                                ][
                                    ()
                                ]  # rotation around G2
                            if rad_range is not None:
                                level_1_subg_3.create_dataset(
                                    "rad_range", dtype="f", data=rad_range
                                )
                            if azim_range is not None:
                                level_1_subg_3.create_dataset(
                                    "azim_range", dtype="f", data=azim_range
                                )
                            if (
                                gon1 is not None
                                or gon2 is not None
                                or gon3 is not None
                                or chiGon1 is not None
                                or omegaGon2 is not None
                                or phiGon3 is not None
                            ):
                                level_1_subg_4 = level_1_subg_1.create_group("Position")
                            print(
                                "### The matrix image is a "
                                + str(numpy.ndim(image))
                                + "D matrix"
                            )
                            flog.write(
                                "### The matrix image is a "
                                + str(numpy.ndim(image))
                                + "D matrix \n"
                            )
                            if (
                                numpy.ndim(image) == 2
                            ):  # do a test on the dimension of the image matrix
                                rslt_matrix_cts = numpy.zeros(
                                    (int(npt_rad), int(npt_azim) + 1), float
                                )
                                rslt_matrix_chi = numpy.zeros((int(npt_azim)), float)
                                rslt_matrix_tth = numpy.zeros((int(npt_rad)), float)
                                if im_mask is None:
                                    mask_mat = None
                                    print("### No mask was used for the integration")
                                    flog.write(
                                        "### No mask was used for the integration \n"
                                    )
                                else:
                                    mask_mat = fabio.open(im_mask).data
                                    print(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask \n"
                                    )
                                if im_dark is None:
                                    dark_mat = None
                                    print("### No dark was used for the integration")
                                    flog.write(
                                        "### No dark was used for the integration \n"
                                    )
                                else:
                                    dark_mat = fabio.open(im_dark).data
                                    print(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark \n"
                                    )
                                if imFlat is None:
                                    flat_mat = None
                                    print(
                                        "### No flat field image was used for the integration"
                                    )
                                    flog.write(
                                        "### No flat field image was used for the integration \n"
                                    )
                                else:
                                    flat_mat = fabio.open(imFlat).data
                                    print(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image \n"
                                    )
                                print(
                                    "### Integration started with a"
                                    + " "
                                    + npt_rad
                                    + " "
                                    + "steps in tth and with a"
                                    + " "
                                    + npt_azim
                                    + " "
                                    + "sectors"
                                )
                                flog.write(
                                    "### Integration started with a"
                                    + " "
                                    + npt_rad
                                    + " "
                                    + "steps in tth and with a"
                                    + " "
                                    + npt_azim
                                    + " "
                                    + "sectors \n"
                                )
                                cts, tth, chi = ai.integrate2d(
                                    image,
                                    int(npt_rad),
                                    int(npt_azim),
                                    correctSolidAngle=True,
                                    error_model=errorModel,
                                    radial_range=rad_range,
                                    azimuth_range=azim_range,
                                    mask=mask_mat,
                                    polarization_factor=0.95,
                                    dark=dark_mat,
                                    flat=flat_mat,
                                    method="splitpixel",
                                    unit=x_unit,
                                )
                                print("### Integration copmleted")
                                flog.write("### Integration copmleted \n")
                                # print(numpy.shape(tth))
                                # print(numpy.shape(cts))
                                # print(numpy.shape(rslt_matrix_cts))
                                rslt_matrix_cts[:, 0] = tth[
                                    :
                                ]  # putting the x axis variables (2tth or q ...) in the first column of the result matrix
                                rslt_matrix_chi[:] = chi[
                                    :
                                ]  # putting the chi angles (around the DS ring) in the chi dataset matrix
                                rslt_matrix_tth[:] = tth[
                                    :
                                ]  # putting the x axis variables (2tth or q ...) in its dataset matrix
                                print("### Saving in h5file started")
                                flog.write("### Saving in h5file started \n")
                                for ichi in range(numpy.shape(chi)[0]):
                                    rslt_matrix_cts[:, ichi + 1] = cts[
                                        ichi, :
                                    ]  # saving the binned intensity as function of chi angles
                                level_1_subg_2 = level_1_subg_1.create_group(
                                    "image_00000"
                                )
                                level_1_subg_2.create_dataset(
                                    "tth_vs_cts", dtype="f", data=rslt_matrix_cts
                                )
                                level_1_subg_2.create_dataset(
                                    "chi", dtype="f", data=rslt_matrix_chi
                                )
                                level_1_subg_2.create_dataset(
                                    "tth", dtype="f", data=rslt_matrix_tth
                                )
                                if gon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S1", dtype="f", data=-g1_axis
                                    )
                                if gon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S2", dtype="f", data=-g2_axis
                                    )
                                if gon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S3", dtype="f", data=-g3_axis
                                    )
                                if chiGon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Chi", dtype="f", data=rg1_chi
                                    )
                                if omegaGon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Omega", dtype="f", data=rg2_omega
                                    )
                                if phiGon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "phi", dtype="f", data=rg3_phi
                                    )
                                print("### Saving in h5file completed")
                                flog.write("### Saving in h5file completed \n")
                            else:  # is executed if the matrix dimension is not 2 (meaninig if it contains more than one image)
                                if im_mask is None:
                                    mask_mat = None
                                    print("### No mask was used for the integration")
                                    flog.write(
                                        "### No mask was used for the integration \n"
                                    )
                                else:
                                    mask_mat = fabio.open(im_mask).data
                                    print(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask \n"
                                    )
                                if im_dark is None:
                                    dark_mat = None
                                    print("### No dark was used for the integration")
                                    flog.write(
                                        "### No dark was used for the integration \n"
                                    )
                                else:
                                    dark_mat = fabio.open(im_dark).data
                                    print(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark \n"
                                    )
                                if imFlat is None:
                                    flat_mat = None
                                    print(
                                        "### No flat field image was used for the integration"
                                    )
                                    flog.write(
                                        "### No flat field image was used for the integration \n"
                                    )
                                else:
                                    flat_mat = fabio.open(imFlat).data
                                    print(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image \n"
                                    )
                                for i in range(numpy.shape(image)[0]):
                                    print(
                                        "*#*#*#*#* Processing of image_"
                                        + "%.0f" % i
                                        + "*#*#*#*#*"
                                    )
                                    flog.write(
                                        "*#*#*#*#* Processing of image_"
                                        + "%.0f" % i
                                        + "*#*#*#*#* \n"
                                    )
                                    rslt_matrix_cts = numpy.zeros(
                                        (int(npt_rad), int(npt_azim) + 1), float
                                    )
                                    rslt_matrix_chi = numpy.zeros(
                                        (int(npt_azim)), float
                                    )
                                    rslt_matrix_tth = numpy.zeros((int(npt_rad)), float)
                                    # mask_mat = numpy.zeros((numpy.shape(image)[0], numpy.shape(image)[1]), float)
                                    # print('### Do not worry I am creating the mask')
                                    # for ix in range(numpy.shape(image)[1]):
                                    # 	for iy in range(numpy.shape(image)[2]):
                                    # 		if (image[i,ix,iy] < numpy.float(int_min) or image[i,ix,iy] > numpy.float(int_max)):
                                    # 			mask_mat[ix,iy] = 1
                                    # print('### Hola, I finished creating the mask')
                                    print(
                                        "### Integration started with a"
                                        + " "
                                        + npt_rad
                                        + " "
                                        + "steps in tth and with a"
                                        + " "
                                        + npt_azim
                                        + " "
                                        + "sectors"
                                    )
                                    flog.write(
                                        "### Integration started with a"
                                        + " "
                                        + npt_rad
                                        + " "
                                        + "steps in tth and with a"
                                        + " "
                                        + npt_azim
                                        + " "
                                        + "sectors \n"
                                    )
                                    cts, tth, chi = ai.integrate2d(
                                        image[i],
                                        int(npt_rad),
                                        int(npt_azim),
                                        correctSolidAngle=True,
                                        error_model=errorModel,
                                        radial_range=rad_range,
                                        azimuth_range=azim_range,
                                        mask=mask_mat,
                                        polarization_factor=0.95,
                                        dark=dark_mat,
                                        flat=flat_mat,
                                        method="splitpixel",
                                        unit=x_unit,
                                    )
                                    print("### Integration copmleted")
                                    flog.write("### Integration copmleted \n")
                                    rslt_matrix_cts[:, 0] = tth[:]
                                    rslt_matrix_chi[:] = chi[:]
                                    rslt_matrix_tth[:] = tth[:]
                                    print("### Saving in h5file started")
                                    flog.write("### Saving in h5file started \n")
                                    for ichi in range(numpy.shape(chi)[0]):
                                        rslt_matrix_cts[:, ichi + 1] = cts[ichi, :]
                                    if i < 10:
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_0000" + "%.0f" % i
                                        )
                                    if (i > 9) and (i < 100):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_000" + "%.0f" % i
                                        )
                                    if (i > 99) and (i < 1000):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_00" + "%.0f" % i
                                        )
                                    if (i > 999) and (i < 10000):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_0" + "%.0f" % i
                                        )
                                    if (i > 9999) and (i < 100000):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_" + "%.0f" % i
                                        )
                                    level_1_subg_2.create_dataset(
                                        "tth_vs_cts", dtype="f", data=rslt_matrix_cts
                                    )
                                    level_1_subg_2.create_dataset(
                                        "chi", dtype="f", data=rslt_matrix_chi
                                    )
                                    level_1_subg_2.create_dataset(
                                        "tth", dtype="f", data=rslt_matrix_tth
                                    )
                                    print("### Saving in h5file completed")
                                    flog.write("### Saving in h5file completed \n")
                                if gon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S1", dtype="f", data=-g1_axis
                                    )
                                if gon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S2", dtype="f", data=-g2_axis
                                    )
                                if gon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S3", dtype="f", data=-g3_axis
                                    )
                                if chiGon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Chi", dtype="f", data=rg1_chi
                                    )
                                if omegaGon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Omega", dtype="f", data=rg2_omega
                                    )
                                if phiGon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "phi", dtype="f", data=rg3_phi
                                    )
                        else:
                            print("### No image was saved in this scan")
                            flog.write("### No image was saved in this scan \n")
                else:  # This loop will be executed if scan='name of scan' and numScan='(start,end)'
                    for i in range(
                        numScan[0], numScan[1] + 1
                    ):  # Iteration on the scans
                        if (
                            "measurement" in r_h5file[scan + str(i) + ".1"]
                            and detector_name
                            in r_h5file[scan + str(i) + ".1/measurement"]
                        ):
                            level_1 = fh5_save.create_group(
                                scan + str(i) + ".1"
                            )  ### Create the group on which will be putted the data of integration
                            level_1_subg_1 = level_1.create_group("raw_integration_2D")
                            level_1_subg_3 = level_1_subg_1.create_group(
                                "Integration_parameter"
                            )
                            level_1_subg_3.create_dataset(
                                "npt_azim", dtype="f", data=int(npt_azim)
                            )
                            level_1_subg_3.create_dataset(
                                "npt_rad", dtype="f", data=int(npt_rad)
                            )
                            image_nc = r_h5file[
                                "/" + scan + str(i) + ".1/measurement/" + detector_name
                            ]  ### Read the image or images
                            image = numpy.float64(
                                image_nc
                            )  # convert the image matrix from int32 to float
                            if gon1 is not None:
                                g1_axis = r_h5file[
                                    "/"
                                    + scan
                                    + str(i)
                                    + ".1"
                                    + "/instrument/positioners/"
                                    + gon1
                                ][
                                    ()
                                ]  # position of the sample along the beam
                            if gon2 is not None:
                                g2_axis = r_h5file[
                                    "/"
                                    + scan
                                    + str(i)
                                    + ".1"
                                    + "/instrument/positioners/"
                                    + gon2
                                ][
                                    ()
                                ]  # position of the samlple in the dir G2
                            if gon3 is not None:
                                g3_axis = r_h5file[
                                    "/"
                                    + scan
                                    + str(i)
                                    + ".1"
                                    + "/instrument/positioners/"
                                    + gon3
                                ][
                                    ()
                                ]  # position of the samlple in the dir G3 (perpendicular to the surface of the sample)
                            if chiGon1 is not None:
                                rg1_chi = r_h5file[
                                    "/"
                                    + scan
                                    + str(i)
                                    + ".1"
                                    + "/instrument/positioners/"
                                    + chiGon1
                                ][
                                    ()
                                ]  # rotation around G1
                            if omegaGon2 is not None:
                                rg2_omega = r_h5file[
                                    "/"
                                    + scan
                                    + str(i)
                                    + ".1"
                                    + "/instrument/positioners/"
                                    + omegaGon2
                                ][
                                    ()
                                ]  # rotation around G2
                            if phiGon3 is not None:
                                rg3_phi = r_h5file[
                                    "/"
                                    + scan
                                    + str(i)
                                    + ".1"
                                    + "/instrument/positioners/"
                                    + phiGon3
                                ][
                                    ()
                                ]  # rotation around G2
                            if rad_range is not None:
                                level_1_subg_3.create_dataset(
                                    "rad_range", dtype="f", data=rad_range
                                )
                            if azim_range is not None:
                                level_1_subg_3.create_dataset(
                                    "azim_range", dtype="f", data=azim_range
                                )
                            if (
                                gon1 is not None
                                or gon2 is not None
                                or gon3 is not None
                                or chiGon1 is not None
                                or omegaGon2 is not None
                                or phiGon3 is not None
                            ):
                                level_1_subg_4 = level_1_subg_1.create_group("Position")
                            print(
                                "### The matrix image is a "
                                + str(numpy.ndim(image))
                                + "D matrix"
                            )
                            flog.write(
                                "### The matrix image is a "
                                + str(numpy.ndim(image))
                                + "D matrix \n"
                            )
                            if (
                                numpy.ndim(image) == 2
                            ):  # do a test on the dimension of the image matrix
                                rslt_matrix_cts = numpy.zeros(
                                    (int(npt_rad), int(npt_azim) + 1), float
                                )
                                rslt_matrix_chi = numpy.zeros((int(npt_azim)), float)
                                rslt_matrix_tth = numpy.zeros((int(npt_rad)), float)
                                if im_mask is None:
                                    mask_mat = None
                                    print("### No mask was used for the integration")
                                    flog.write(
                                        "### No mask was used for the integration \n"
                                    )
                                else:
                                    mask_mat = fabio.open(im_mask).data
                                    print(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask \n"
                                    )
                                if im_dark is None:
                                    dark_mat = None
                                    print("### No dark was used for the integration")
                                    flog.write(
                                        "### No dark was used for the integration \n"
                                    )
                                else:
                                    dark_mat = fabio.open(im_dark).data
                                    print(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark \n"
                                    )
                                if imFlat is None:
                                    flat_mat = None
                                    print(
                                        "### No flat field image was used for the integration"
                                    )
                                    flog.write(
                                        "### No flat field image was used for the integration \n"
                                    )
                                else:
                                    flat_mat = fabio.open(imFlat).data
                                    print(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image \n"
                                    )
                                print(
                                    "### Integration started with a"
                                    + " "
                                    + str(npt_rad)
                                    + " "
                                    + "steps in tth and with a"
                                    + " "
                                    + str(npt_azim)
                                    + " "
                                    + "sectors"
                                )
                                flog.write(
                                    "### Integration started with a"
                                    + " "
                                    + str(npt_rad)
                                    + " "
                                    + "steps in tth and with a"
                                    + " "
                                    + str(npt_azim)
                                    + " "
                                    + "sectors \n"
                                )
                                cts, tth, chi = ai.integrate2d(
                                    image,
                                    int(npt_rad),
                                    int(npt_azim),
                                    correctSolidAngle=True,
                                    error_model=errorModel,
                                    radial_range=rad_range,
                                    azimuth_range=azim_range,
                                    mask=mask_mat,
                                    polarization_factor=0.95,
                                    dark=dark_mat,
                                    flat=flat_mat,
                                    method="splitpixel",
                                    unit=x_unit,
                                )
                                print("### Integration copmleted")
                                flog.write("### Integration copmleted \n")
                                # print(numpy.shape(tth))
                                # print(numpy.shape(cts))
                                # print(numpy.shape(rslt_matrix_cts))
                                rslt_matrix_cts[:, 0] = tth[
                                    :
                                ]  # putting the x axis variables (2tth or q ...) in the first column of the result matrix
                                rslt_matrix_chi[:] = chi[
                                    :
                                ]  # putting the chi angles (around the DS ring) in the chi dataset matrix
                                rslt_matrix_tth[:] = tth[
                                    :
                                ]  # putting the x axis variables (2tth or q ...) in its dataset matrix
                                print("### Saving in h5file started")
                                flog.write("### Saving in h5file started \n")
                                for ichi in range(numpy.shape(chi)[0]):
                                    rslt_matrix_cts[:, ichi + 1] = cts[
                                        ichi, :
                                    ]  # saving the binned intensity as function of chi angles
                                level_1_subg_2 = level_1_subg_1.create_group(
                                    "image_00000"
                                )
                                level_1_subg_2.create_dataset(
                                    "tth_vs_cts", dtype="f", data=rslt_matrix_cts
                                )
                                level_1_subg_2.create_dataset(
                                    "chi", dtype="f", data=rslt_matrix_chi
                                )
                                level_1_subg_2.create_dataset(
                                    "tth", dtype="f", data=rslt_matrix_tth
                                )
                                if gon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S1", dtype="f", data=-g1_axis
                                    )
                                if gon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S2", dtype="f", data=-g2_axis
                                    )
                                if gon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S3", dtype="f", data=-g3_axis
                                    )
                                if chiGon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Chi", dtype="f", data=rg1_chi
                                    )
                                if omegaGon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Omega", dtype="f", data=rg2_omega
                                    )
                                if phiGon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "phi", dtype="f", data=rg3_phi
                                    )
                                print("### Saving in h5file completed")
                                flog.write("### Saving in h5file completed \n")
                            else:  # is executed if the matrix dimension is not 2 (meaninig if it contains more than one image)
                                if im_mask is None:
                                    mask_mat = None
                                    print("### No mask was used for the integration")
                                    flog.write(
                                        "### No mask was used for the integration \n"
                                    )
                                else:
                                    mask_mat = fabio.open(im_mask).data
                                    print(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_mask
                                        + " "
                                        + "was used as mask \n"
                                    )
                                if im_dark is None:
                                    dark_mat = None
                                    print("### No dark was used for the integration")
                                    flog.write(
                                        "### No dark was used for the integration \n"
                                    )
                                else:
                                    dark_mat = fabio.open(im_dark).data
                                    print(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + im_dark
                                        + " "
                                        + "was used as dark \n"
                                    )
                                if imFlat is None:
                                    flat_mat = None
                                    print(
                                        "### No flat field image was used for the integration"
                                    )
                                    flog.write(
                                        "### No flat field image was used for the integration \n"
                                    )
                                else:
                                    flat_mat = fabio.open(imFlat).data
                                    print(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image"
                                    )
                                    flog.write(
                                        "### The image: "
                                        + imFlat
                                        + " "
                                        + "was used as flat field image \n"
                                    )
                                for i in range(numpy.shape(image)[0]):
                                    print(
                                        "*#*#*#*#* Processing of image_"
                                        + "%.0f" % i
                                        + "*#*#*#*#*"
                                    )
                                    flog.write(
                                        "*#*#*#*#* Processing of image_"
                                        + "%.0f" % i
                                        + "*#*#*#*#* \n"
                                    )
                                    rslt_matrix_cts = numpy.zeros(
                                        (int(npt_rad), int(npt_azim) + 1), float
                                    )
                                    rslt_matrix_chi = numpy.zeros(
                                        (int(npt_azim)), float
                                    )
                                    rslt_matrix_tth = numpy.zeros((int(npt_rad)), float)
                                    # mask_mat = numpy.zeros((numpy.shape(image)[0], numpy.shape(image)[1]), float)
                                    # print('### Do not worry I am creating the mask')
                                    # for ix in range(numpy.shape(image)[1]):
                                    # 	for iy in range(numpy.shape(image)[2]):
                                    # 		if (image[i,ix,iy] < numpy.float(int_min) or image[i,ix,iy] > numpy.float(int_max)):
                                    # 			mask_mat[ix,iy] = 1
                                    # print('### Hola, I finished creating the mask')
                                    print(
                                        "### Integration started with a"
                                        + " "
                                        + str(npt_rad)
                                        + " "
                                        + "steps in tth and with a"
                                        + " "
                                        + str(npt_azim)
                                        + " "
                                        + "sectors"
                                    )
                                    flog.write(
                                        "### Integration started with a"
                                        + " "
                                        + str(npt_rad)
                                        + " "
                                        + "steps in tth and with a"
                                        + " "
                                        + str(npt_azim)
                                        + " "
                                        + "sectors \n"
                                    )
                                    cts, tth, chi = ai.integrate2d(
                                        image[i],
                                        int(npt_rad),
                                        int(npt_azim),
                                        correctSolidAngle=True,
                                        error_model=errorModel,
                                        radial_range=rad_range,
                                        azimuth_range=azim_range,
                                        mask=mask_mat,
                                        polarization_factor=0.95,
                                        dark=dark_mat,
                                        flat=flat_mat,
                                        method="splitpixel",
                                        unit=x_unit,
                                    )
                                    print("### Integration copmleted")
                                    flog.write("### Integration copmleted \n")
                                    rslt_matrix_cts[:, 0] = tth[:]
                                    rslt_matrix_chi[:] = chi[:]
                                    rslt_matrix_tth[:] = tth[:]
                                    print("### Saving in h5file started")
                                    flog.write("### Saving in h5file started \n")
                                    for ichi in range(numpy.shape(chi)[0]):
                                        rslt_matrix_cts[:, ichi + 1] = cts[ichi, :]
                                    if i < 10:
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_0000" + "%.0f" % i
                                        )
                                    if (i > 9) and (i < 100):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_000" + "%.0f" % i
                                        )
                                    if (i > 99) and (i < 1000):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_00" + "%.0f" % i
                                        )
                                    if (i > 999) and (i < 10000):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_0" + "%.0f" % i
                                        )
                                    if (i > 9999) and (i < 100000):
                                        level_1_subg_2 = level_1_subg_1.create_group(
                                            "image" + "_" + "%.0f" % i
                                        )
                                    level_1_subg_2.create_dataset(
                                        "tth_vs_cts", dtype="f", data=rslt_matrix_cts
                                    )
                                    level_1_subg_2.create_dataset(
                                        "chi", dtype="f", data=rslt_matrix_chi
                                    )
                                    level_1_subg_2.create_dataset(
                                        "tth", dtype="f", data=rslt_matrix_tth
                                    )
                                    print("### Saving in h5file completed")
                                    flog.write("### Saving in h5file completed \n")
                                if gon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S1", dtype="f", data=-g1_axis
                                    )
                                if gon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S2", dtype="f", data=-g2_axis
                                    )
                                if gon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "S3", dtype="f", data=-g3_axis
                                    )
                                if chiGon1 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Chi", dtype="f", data=rg1_chi
                                    )
                                if omegaGon2 is not None:
                                    level_1_subg_4.create_dataset(
                                        "Omega", dtype="f", data=rg2_omega
                                    )
                                if phiGon3 is not None:
                                    level_1_subg_4.create_dataset(
                                        "phi", dtype="f", data=rg3_phi
                                    )
                        else:
                            print("### No image was saved in this scan")
                            flog.write("### No image was saved in this scan \n")
                print("### Hope to see you again")
                flog.write("### Hope to see you again \n")
