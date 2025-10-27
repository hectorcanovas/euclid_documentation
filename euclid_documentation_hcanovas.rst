.. _astroquery.esa.euclid:

********************************************
ESA EUCLID Archive (`astroquery.esa.euclid`)
********************************************


**The Euclid mission**

`Euclid <https://www.cosmos.esa.int/web/euclid>`_ is a Medium class ESA mission to map the geometry of the dark Universe.
The mission investigates the distance-redshift relationship and the evolution of cosmic structures. The space telescope is
creating a detailed map of the large-scale structure of the Universe across space and time by observing billions of galaxies
out to 10 billion light-years, across more than a third of the sky. It achieves this by measuring shapes and redshifts of
galaxies and clusters of galaxies out to redshifts ~2, or equivalently to a look-back time of 10 billion years. It therefore
explores how the Universe has expanded and how structure has formed over cosmic history, revealing more about the role of
gravity and the nature of dark energy and dark matter.

Please take note of our `guide <https://www.cosmos.esa.int/web/euclid/data-credits-acknowledgements>`_ on how
to acknowledge and cite Euclid data if you use public Euclid data in your paper.

**hcanovas: remove this paragraph. WHY: it provides too much detail about the mission (compared to other Astroquery pages).
Proposal: add it to the Euclid Archive Help pages (location to be decided):**
The `Euclid Survey <https://www.euclid-ec.org/public/data/surveys/>`_ is done in a 'step-and-stare' mode, where the telescope
points to a position on the sky and then imaging and spectroscopic measurements are performed on an area of ~0.48 deg\ :sup:`2`
around this position. The telescope consists of two cameras, the visible instrument (VIS) and the Near Infrared Spectrometer
and Photometer (NISP) instrument that observe simultaneously using a light splitting dichroic.
For the survey standard operating mode, the telescope undertakes a 4-point dither pattern. At each position VIS and NISP
each take a 560s exposure, consisting of a direct visible image and a red grism exposure. This is followed by further
NISP exposures in the Y, J, and H band filters (87 seconds each). The telescope is then dithered, and the sequence is
repeated starting with a different grism position angle. There are actually two operational grisms oriented 180 degrees
from each other. Each grism will be used twice in this sequence, but with slight angular offsets (+/- 4 degrees),
effectively creating the four different grism angles
(`Scaramella et al. 2022, A&A 662, A112 <https://ui.adsabs.harvard.edu/abs/2022A%26A...662A.112E/abstract>`_).
This standard four-dithers operating mode sequence is called a single observation and all the individual exposures
associated with each observation are organized by Observation ID in the archive. The
`Science Ground Segment <https://www.euclid-ec.org/public/data/ground-segment/>`_ also processes all of its
imaging into merged mosaics, which can contain multiple different observations. All products associated with these
mosaics are organized by Tile ID in the archive.

**hcanovas: remove the following four paragraphs (and example). WHY: it is not intended for public users. None of the astroquery pages provides info like this.
Proposal: Add it to the EC Help content & EC datalabs:** The Euclid Science Archive has several environments
serving different purposes for the `Euclid Consortium <https://www.euclid-ec.org/>`_ members.

1. The OTF ("on-the-fly") environment of the Euclid science archive, first started at the start of science operation exposed data as processed by
the SGS (Science Ground Segment) soon after acquisition to provide an access as soon as possible. In this environment
the data will not be reprocessed and the processing is therefore heterogeneous.

2. The REG (for non-regression testing) environment of the Euclid science archive, where a large area in the sky is processed with the same version for
all data products. The first campaign was run in September 2024, for area of about 500 square degrees (~1000
observations), the next campaign shall be run in March-April 2025.

3. The IDR (Internal Data Release) environment of the Euclid science archive holds the data that will then become public. The first release Q1
opened on the 6th of November 2024, with a first pass on the three Euclid deep fields (EDFN, EDFS and EDFF) as well as
observations on the Lynds Dark Nebula LDN1641.

4. The PDR (Public Data Release) environment of the Euclid science archive holds the public data. Euclid Q1_ data was publicly released on March 19,
2025. The main component of the Q1_ data contains Level 2 data of a single visit (at the depth of the Euclid Wide
Survey) over the Euclid Deep Fields (EDFs): 20 deg\ :sup:`2` of the EDF North, 10 deg\ :sup:`2` of EDF Fornax, and
23 deg\ :sup:`2` of the EDF South. The deep fields will be visited multiple times during the mission.


By default, the object *Euclid*

  >>> from astroquery.esa.euclid import Euclid

makes use of the *PDR* environment. In order to make use of a different one, it is necessary to instantiate the class EuclidClass

  >>> from astroquery.esa.euclid import EuclidClass
  >>> euclid = EuclidClass(environment='IDR')

The parameter *environment* is limited to *IDR*, *OTF*, *PDR* or *REG*.


**Astroquery.esa.euclid**

This Python module provides an Astroquery API to access to the metadata and datasets provided by the 
`European Space Agency Euclid Archive <https://eas.esac.esa.int/sas/>`_ using a `TAP+ <https://astroquery.readthedocs.io/en/latest/utils/tap.html>`_ REST service.
`TAP+ <https://astroquery.readthedocs.io/en/latest/utils/tap.html>`_ is an extension of Table Access Protocol (TAP_)
specified by the International Virtual Observatory Alliance (IVOA_) that incorporates dedicated user space capabilities (see Sect. 2 below). 
The TAP_ query language is Astronomical Data Query Language (ADQL_). TAP_ provides two operation modes: 

* Synchronous: the server response is generated as soon as the request is received.

* Asynchronous: the server starts a job that will execute the request. The first response to the request is a link with information about the job status.


On top of that, this package provides two access modes: 

* Public (default): The results generated by the anonymous ADQL_ queries are public, and deleted from the Archive 72 hours.

* Authenticated: The ADQL_ queries and their outcomes remain in the user space until the user deletes them. In addition, authenticated users benefit from dedicated functionalities (see Sect 2 below).

There are limitations to the execution time and total output size that depend on the combination of operation and access modes - see the Gaia Archive FAQ:
`Why does my query time out after 90 minutes? Why is my query limited to 3 million rows? <https://www.cosmos.esa.int/web/gaia/faqs#account-limits-2020>`_ for details.


To reduce the examples verbosity (as well as complexity), the code examples output has been trimmed
and only the most relevant output lines are displayed. Whenever possible, the documentation points
to the `Astroquery.Gaia package <https://astroquery.readthedocs.io/en/latest/gaia/gaia.html>`_ that
shares a similar architecture and methods with this module.



**Euclid data and data access**

Euclid Q1_ contains different types of data, like catalogues (data tables), images, and spectra. For details, please refer to the `Data products in the science archive <https://s2e2.cosmos.esa.int/www/euclid_iscience/Data_products_in_the_science_archive.html>`_ in the `Euclid Archive Help <https://s2e2.cosmos.esa.int/www/euclid_iscience/Public_User_Guide.html>`_ , as well as the Q1 Data Product Definition Document (DPDD_). 

This Astroquery package is mostly geared to query and retrieve the data stored in the catalogues, but it also includes dedicated methods to retrieve the images and spectra (both stored as large FITS files). The latter are served via the DataLink_ IVOA_ protocol - see Sect. 3 below). It is also possible to directly access to these products (without having to retrieve them) using the "Euclid Q1" datalab that is publicly available in the ESA Datalabs_ e-science platform. Users aiming to analyse large Euclid datasets are encouraged to use this platform.


Table of contents:

.. contents::
   :local:
   :depth: 3

========
Examples
========

It is recommended checking the status of Euclid TAP_ before executing this module. To do this:

.. almost all code examples require remote-data access, thus only using this
   one at the first example
.. doctest-remote-data-all::

  >>> from astroquery.esa.euclid import Euclid
  >>> Euclid.get_status_messages()


1. Public access
---------------------------

1.1. Metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Table and column metadata are specified by the IVOA_ TAP_ recommendation. To load only table names:

  >>> tables = Euclid.load_tables(only_names=True, include_shared_tables=True)
  >>> print(f'* Found {len(tables)} tables')
  >>> print(*(table.name for table in tables), sep="\n")  # doctest: +IGNORE_OUTPUT
  ivoa.obscore
  public.dual
  sedm.raw_detector
  sedm.raw_frame
  sedm.raw_quadrant
  ...


To load all tables metadata:

  >>> tables = Euclid.load_tables()
  >>> print(tables[0])
  TAP Table name: ivoa.obscore
  Description: None
  Size (bytes): 0
  Num. columns: 34


To load only one table and inspect its columns:

  >>> raw_detector_table = Euclid.load_table('sedm.raw_detector')
  >>> print(raw_detector_table)
  TAP Table name: sedm.raw_detector
  Description: None
  Size (bytes): 0
  Num. columns: 12
  >>> print(*(column.name for column in raw_detector_table.columns), sep="\n")  # doctest: +IGNORE_OUTPUT
  crpix1
  crpix2
  crval1
  ...


1.2. Cone search
^^^^^^^^^^^^^^^^

The cone_search_ method implements one of the most popular use cases when connecting to an astronomy archive: retrieving data around a projected circular region in a given sky location from a given catalogue.
The example below shows how to launch a 0.5 degrees radius cone search around `NGC 6505 <https://simbad.cds.unistra.fr/simbad/sim-id?Ident=NGC+6505>`_. By default, this method targets
the "mer_catalogue" and its outcome is restricted to 50 rows. This limitation can be removed by setting the ROW_LIMIT attribute to "-1" (see the next example below).

  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> coord  = SkyCoord("17h51m07.4s +65d31m50.8s", frame='icrs')
  >>> radius = u.Quantity(0.5, u.deg)
  >>> job    = Euclid.cone_search(coordinate=coord, radius=radius, columns="*", async_job=True)
  >>> res    = job.get_results()
  >>> print(f"Found {len(cone_results)} results")
  basic_download_data_oid to_be_published      object_id       right_ascension   ...       gaia_id        gaia_match_quality           dist         
  ----------------------- --------------- ------------------- ------------------ ... ------------------- -------------------- ----------------------
                      281               1 2677813028655307424 267.78130284070573 ...                  --                   -- 0.00019012520229516453
                      281               1 2677926210655368830  267.7926210570132 ... 1441085261522268928  0.05010449141263962   0.007807472554918013
                      281               1 2677747417655202562  267.7747417649051 ... 1441085055363835648 0.002425010548904538   0.010827275337244202


The example below shows how to 1) remove the row limitation, and 2) target a different table. It also shows that the cone_search method accepts target names of coordinates, provided
that the name is recognised by the Simbad, VizieR, or NED services.

  >>> Euclid.ROW_LIMIT = -1   # Set this attribute to -1 to retrieve the full cone search output.
  >>> job              = Euclid.cone_search(coordinate='NGC 6505', radius=radius, table_name="sedm.mosaic_product", ra_column_name="ra", dec_column_name="dec", columns="*", async_job=True)
  >>> results          = job.get_results()
  >>> print(results[0:3])
  category             checksum                  creation_date      crpix1 crpix2 ... tile_index to_be_published zero_point zero_point_error        dist       
  -------- -------------------------------- ----------------------- ------ ------ ... ---------- --------------- ---------- ---------------- ------------------
  SCIENCE 528dcb14904e7501fca6f2cebf112f38 2024-10-26T14:01:21.038 9600.0 9720.0 ...  102158889               1       30.0              0.1 0.1689677160687657
  SCIENCE 77a35773063c6a088e92b294df817e7e 2024-10-26T13:50:13.676 9600.0 9720.0 ...  102158889               1       30.0              0.1 0.1689677160687657
  SCIENCE 758ae0c3d04d21544facd5d241ae3a07 2024-10-26T13:37:09.628 9600.0 9720.0 ...  102158889               1       29.8              0.1 0.1689677160687657
  ...

**Notes:**

* Once the table_name, and/or ra_column_name, and/or dec_column_name arguments are set, the default values are erased - this is a known issue.

* Users are encouraged to use the cone_search_ instead of the query_object_ method. The latter makes use of the ADQL_ BOX function that is deprecated and can yield misleading results due to geometric projection effects. 





1.3. Synchronous query
^^^^^^^^^^^^^^^^^^^^^^

This is the recommended access mode for queries that do not require excessive computation time and/or generate tables with less than 2,000 rows - for details please see the Gaia Archive FAQ:
`Why does my query time out after 90 minutes? Why is my query limited to 3 million rows? <https://www.cosmos.esa.int/web/gaia/faqs#account-limits-2020>`_.
The example below shows how to extract a subset of three sources with ellipticity larger than zero from the "mer_catalogue":


  >>> query = f"SELECT TOP 3 object_id, right_ascension, declination, segmentation_area, ellipticity, kron_radius FROM {mer_cat_name} WHERE ellipticity > 0"
  >>> job   = Euclid.launch_job(query)
  >>> res   = job.get_results()
  >>> print(res)
       object_id       right_ascension      declination    segmentation_area     ellipticity        kron_radius    
  ------------------- ------------------ ----------------- ----------------- ------------------- ------------------
  2744182404684288509 274.41824043555044 68.42885096091729                98 0.34178537130355835 22.018617630004883
  2744679115684290125  274.4679115369266  68.4290125712924                94 0.36368849873542786  19.16196632385254
  2744820013684293317 274.48200131993156 68.42933179142987                47 0.13406670093536377 13.094022750854492


The launch_job_ method returns a *Job* object. Its results can be extracted using the "get_results()" method, that generates an `Astropy table <https://docs.astropy.org/en/stable/table/index.html>`_ object.
The job status can be inspected by typing:

  >>> print(job)

Note that deleting the "TOP 3" string in the query above will return a table with 2,000 rows (sources).


1.4. Asynchronous query
^^^^^^^^^^^^^^^^^^^^^^^

This is the recommended mode for queries that are expected to output more than 2,000 rows and that require substantial execution time
(noting that all the queries time out after 7200 seconds). The query results are stored in the Archive (during 72 hours for anonymous users, and on the user area until the user deletes them for registered users). The example below generates a cone search combined with a constraint applied to the ellipticity, and is similar to the first ADQL_ query example listed in the `Euclid Archive <https://eas.esac.esa.int/sas/>`_ (see its "Search/ADQL FORM" subtab). For more ADQL_ examples please have a look at the Gaia Archive Help content (in particular, the `writting queries <https://www.cosmos.esa.int/web/gaia-users/archive/writing-queries>`_ section). 


  >>> query = "SELECT right_ascension, declination, object_id, vis_det, det_quality_flag, flux_detection_total, flux_vis_sersic, segmentation_area, kron_radius, DISTANCE(267.78, 65.53, right_ascension, declination) AS dist FROM mer_catalogue WHERE DISTANCE(267.78, 65.53, right_ascension, declination) < 0.1 AND ellipticity > 0"
  >>> job  = Euclid.launch_job_async(query, verbose=False)
  >>> print(job_async)
  >>> res  = job.get_results()
  >>> print(res)  
   right_ascension      declination         object_id      vis_det det_quality_flag ...   flux_vis_sersic   segmentation_area    kron_radius             dist       
  ------------------ ----------------- ------------------- ------- ---------------- ... ------------------- ----------------- ------------------ -------------------
   267.7502407456637 65.43182123675119 2677502407654318212       1                2 ...  0.5517764091491699               101 28.096778869628906 0.09895246869367581
  267.76847561971346 65.43194661918689 2677684756654319466       1                2 ... 0.12194870412349701                26 12.046527862548828  0.0981699461580744
   267.7698066292422 65.43454030481347 2677698066654345403       1                0 ...  0.4004257917404175                61 15.562020301818848 0.09555336819768735
   ...



1.5 Query on an 'on-the-fly' uploaded table
^^^^^^^^^^^^^^^^^^^^^^^

This feature is present both in the synchronous and asynchronous requests. 'On-the-fly' queries allow you to upload a table stored in VOTable_ format and
perform a query on it in one single command. The uploaded tables are deleted after the query is complete. Alternatively, as a registered user it is possible
to upload a table and store it in the user space (see Sect. 2 below). In the example below, the "my_table.xml" file is uploaded to the Archive and used to
perform a JOIN operation with the mer_catalogue. Note the use of the "tap_upload" in the ADQL_ query.

  >>> upload_resource = 'my_table.xml'
  >>> query           = "SELECT mer.object_id, flux_vis_sersic, fwhm FROM tap_upload.table_test JOIN mer_catalogue AS mer USING (object_id)"
  >>> job             = Euclid.launch_job(query, upload_resource=upload_resource, upload_table_name="table_test")
  >>> res             = job.get_results()
  >>> print(res)
        object_id        flux_vis_sersic          fwhm       
   ------------------- ------------------- ------------------
   2701338214642376775 0.12060821056365967 1.2710362672805786
   2703077159642376093  0.8660399317741394 1.2481366395950317
   2695939228642370900    6.01658296585083 1.2056430578231812



1.6. Getting products
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A popular Euclid use case is "I want to download the Euclid image of my favourite object". In order to do this, one must know the coordinates of the selected object as well as the type of image (e.g. calibrated, mosaic, ...) and instrument that was used to acquire the image (e.g., VIS)
that is going to be downloaded. In the following example we use the get_product_ method to download one VIS calibrated image (the Euclid product: DpdVisCalibratedQuadFrame_) that contains the sky region where NGC 6505 is located). The image is then displayed using Matplotlib. 


**Notes:**
* The Euclid FITS files are large (~7 GB for calibrated VIS images)
* When running this notebook on ESA Datalabs you have access to all the products directly through the data volume and don't need to download them - this method is an alternative for working outside ESA Datalabs or accessing data that is not in a data volume.


.. _TAP: http://www.ivoa.net/documents/TAP/
.. _IVOA: http://www.ivoa.net
.. _ADQL: https://www.ivoa.net/documents/ADQL/20231215/index.html
.. _DataLink: https://www.ivoa.net/documents/DataLink/20231215/index.html
.. _VOTable: https://www.ivoa.net/documents/VOTable/20250116/
.. _Q1: https://www.cosmos.esa.int/web/euclid/euclid-q1-data-release
.. _DPDD: https://euclid.esac.esa.int/dr/q1/dpdd/index.html
.. _REST: https://en.wikipedia.org/wiki/Representational_state_transfer
.. _cone_search: https://astroquery.readthedocs.io/en/latest/api/astroquery.esa.euclid.EuclidClass.html#astroquery.esa.euclid.EuclidClass.cone_search
.. _query_object: https://astroquery.readthedocs.io/en/latest/api/astroquery.esa.euclid.EuclidClass.html#astroquery.esa.euclid.EuclidClass.query_object
.. _launch_job: https://astroquery.readthedocs.io/en/latest/api/astroquery.esa.euclid.EuclidClass.html#astroquery.esa.euclid.EuclidClass.launch_job 
.. _launch_job_async: https://astroquery.readthedocs.io/en/latest/api/astroquery.esa.euclid.EuclidClass.html#astroquery.esa.euclid.EuclidClass.launch_job_async 
.. _get_product: https://astroquery.readthedocs.io/en/latest/api/astroquery.esa.euclid.EuclidClass.html#astroquery.esa.euclid.EuclidClass.get_product
.. _DpdVisCalibratedQuadFrame: https://euclid.esac.esa.int/dr/q1/dpdd/visdpd/dpcards/vis_calibratedquadframe.html
.. _Datalabs: https://datalabs.esa.int/


It is possible to download a product given its file name or product id:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> #makeing a folder for the output files
  >>> import os
  >>> output_folder= 'example_outputs/'
  >>> if not os.path.exists(output_folder):
         os.makedirs(output_folder)
  >>>
  >>> example_file_name = "EUC_MER_BGSUB-MOSAIC-NIR-H_TILE102158889-ED035A_20241024T212936.705156Z_00.00.fits"
  >>> print("Getting file:", example_file_name)
  Getting file: EUC_MER_BGSUB-MOSAIC-NIR-H_TILE102158889-ED035A_20241024T212936.705156Z_00.00.fits
  >>> path = Euclid.get_product(file_name=example_file_name, output_file=output_folder + example_file_name,verbose=True)
  Retrieving data.
  Data request: TAPCLIENT=ASTROQUERY&RELEASE=sedm&FILE_NAME=EUC_MER_BGSUB-MOSAIC-NIR-H_TILE102158889-ED035A_20241024T212936.705156Z_00.00.fits&RETRIEVAL_TYPE=FILE
  ------>https
  host = easidr.esac.esa.int:443
  context = /sas-dd/data
  Content-type = application/x-www-form-urlencoded
  200
  Reading...
  Done.
  >>> #display the downloaded product (since this is a calibrated frame the different detectors are stored as different extensions - we are displaying only one extension)
  >>> from astropy.io import fits
  >>> import matplotlib.pyplot as plt
  >>> from astropy.visualization import astropy_mpl_style, ImageNormalize, PercentileInterval, AsinhStretch, LogStretch
  >>> hdul = fits.open(path[0])
  >>> print(fits.info(path[0]))
  WARNING: File may have been truncated: actual file length (103579232) is smaller than the expected size (1474565760) [astropy.io.fits.file]
  Filename: example_notebook_outputs/EUC_MER_BGSUB-MOSAIC-DES-I_TILE102018211-31E2C9_20241018T143048.358037Z_00.00.fits
  No.    Name      Ver    Type      Cards   Dimensions   Format
    0  PRIMARY       1 PrimaryHDU      48   (19200, 19200)   float32
  None
  >>> image_data = hdul[0].data
  >>>
  >>> plt.figure()
  <Figure size 800x600 with 0 Axes>
  <Figure size 800x600 with 0 Axes>
  >>> plt.imshow(image_data, cmap='gray', origin='lower', norm=ImageNormalize(image_data, interval=PercentileInterval(99.9), stretch=AsinhStretch()))
  >>> colorbar = plt.colorbar()


.. image:: images/EUC_MER_BGSUB-MOSAIC-NIR-H_TILE102158889-ED035A_20241024T212936.705156Z_00.00.png
   :align: center
   :scale: 100%
   :alt: EUC_MER_BGSUB-MOSAIC-NIR-H_TILE102158889-ED035A_20241024T212936.705156Z_00.00.fits


The method downloads the fits file(s) and returns the local path where the product(s) is saved.

To download the products for a given EUCLID observation_id (observations) or tile_index (mosaics):

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> #downloading all products for observation id: 102018211
  >>> from astroquery.esa.euclid import Euclid
  >>> mos_id = 1399
  >>> path = Euclid.get_observation_products(id=mos_id, product_type='mosaic', filter="VIS", output_file=f"{output_folder}/products_{mos_id}.fits", verbose=True)

  For big files the download may require a long time.


1.7. Cutout search
^^^^^^^^^^^^^^^^^^

To download a cutout given its file path, instrument and obs_id, and the cutout region, the method downloads the fits file of the cutout and returns a list containing the local path where the cutout is saved:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> # the map cone_results was previously obtained by the query executed in section 2.1
  >>> from astroquery.esa.euclid import Euclid
  >>> from astropy.coordinates import SkyCoord
  >>> import astropy.units as u
  >>> example_file = cone_results[cone_results['instrument_name'] == 'VIS'][0]
  >>> # getting the arguments from the cone search result table automatically
  >>> file_path=example_file["file_path"] + "/" + example_file["file_name"]
  >>> instrument=example_file["instrument_name"]
  >>> obs_id=example_file["tile_index"]
  >>> radius= 0.2 * u.arcmin
  >>> coord = SkyCoord("17h51m07.4s +65d31m50.8s", frame='icrs')
  >>> output_folder= 'example_outputs/'
  >>> if not os.path.exists(output_folder):
         os.makedirs(output_folder)
  >>> output_file=output_folder + 'cutouts/astroquery_cutout_example.fits'
  >>> saved_cutout_filepath = Euclid.get_cutout(file_path=file_path, instrument=instrument, id=obs_id, coordinate=coord, radius=radius, output_file=output_file)
  >>> print("Cutout saved at", saved_cutout_filepath)
  Cutout saved at ['example_outputs/cutouts/astroquery_cutout_example.fits']
  >>>
  >>> #looking at the cutout we made
  >>> hdul = fits.open(saved_cutout_filepath[0])
  >>> print(fits.info(saved_cutout_filepath[0]))
  Filename: example_notebook_outputs/cutouts/astroquery_cutout_example.fits
  >>> image_data = hdul[0].data
  No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      49   (241, 241)   float32
  None
  >>> plt.imshow(image_data, interpolation='nearest', cmap='gray', origin='lower', norm=ImageNormalize(image_data, interval=PercentileInterval(99.5), stretch=AsinhStretch()))
  >>> plt.colorbar()



.. image:: images/astroquery_cutout_example.png
   :align: center
   :scale: 100%
   :alt: astroquery_cutout_example.fits


Below is the equivalent version but copying arguments manually (for clarity).

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> file_path="EUC_MER_BGSUB-MOSAIC-VIS_TILE101158889-D08FBD_20240113T021028.995617Z_00.00.fits"
  >>> saved_cutout_filepath = Euclid.get_cutout(file_path=file_path, instrument="VIS", id="101158889", coordinate=coord, radius=radius, output_file='example_outputs/test_cutout_example.fits')
  >>> print("Cutout saved at", saved_cutout_filepath)
  Cutout saved at ['example_outputs/cutouts/astroquery_cutout_example.fits']






1.8. Getting product data (only useful for DR1 (I think...)!)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the list of products associated with a given Euclid observation_id or tile_index (for mosaic):

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> product_list_results = Euclid.get_product_list(tile_index="102018211", product_type="DpdMerBksMosaic")
  >>> print("Found", len(product_list_results), "results")
  Found 12 results
  >>> print(product_list_results)
                                      file_name                                      mosaic_product_oid tile_index instrument_name filter_name category second_type     ra       dec   technique
                                        str255                                             int64          int64         str255        str255    str255     str255    float64   float64   str255
  ---------------------------------------------------------------------------------- ------------------ ---------- --------------- ----------- -------- ----------- ---------- ------- ---------
  EUC_MER_BGSUB-MOSAIC-DES-I_TILE102018211-31E2C9_20241018T143048.358037Z_00.00.fits               1399  102018211           DECAM     DECAM_i  SCIENCE         SKY 57.9990741   -51.5     IMAGE
    EUC_MER_BGSUB-MOSAIC-VIS_TILE102018211-ACBD03_20241018T142710.276838Z_00.00.fits               1395  102018211             VIS         VIS  SCIENCE         SKY 57.9990741   -51.5     IMAGE
  EUC_MER_BGSUB-MOSAIC-DES-G_TILE102018211-D9D163_20241018T143010.768685Z_00.00.fits               1394  102018211           DECAM     DECAM_g  SCIENCE         SKY 57.9990741   -51.5     IMAGE
  ...

The method returns a list of products as an `~astropy.table.Table`. It is also possible to search by observation_id, but not by both parameters simultaneously.

It is possible to retrieve LE3 data (scientific data) by observation_id or tile_index (but not by both simultaneously) and/or for different categories, groups and product types. The available values
for these parameters are summarized in section :ref:`appendix`.


.. Skipping authentication requiring examples
.. doctest-skip::

  >>> le3_product_list = Euclid.get_scientific_product_list(tile_index=22)
  >>> print("Found", len(le3_product_list), "results")
  Found 3 results
  >>> print(le3_product_list)
  basic_download_data_oid  product_type                            product_id                          observation_id_list tile_index_list patch_id_list filter_name
  ----------------------- -------------- ------------------------------------------------------------- ------------------- --------------- ------------- -----------
                    47191 DpdLE3clCLTile       PPO_REGREPROC1_R2_CLTEST_R0_CLTILING_R5-output_tiles-27                  {}            {22}            {}
                    47132 DpdLE3clCLTile PPO_REGREPROC1_R2_CLTEST_R0_CLTILINGPOLYHR_R2-output_tiles-27                  {}            {22}            {}
                    47233 DpdLE3clCLTile       PPO_REGREPROC1_R2_CLTEST_R0_CLTILING_R6-output_tiles-27                  {}            {22}            {}


In the following example, for the Clusters of Galaxies category, and the group GrpCatalog, we retrieve all the DET-CL AMICO auxiliary Data Product products (DpdLE3clAmicoAux):

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> results = euclid.get_scientific_product_list(category='Clusters of Galaxies', group='GrpCatalog', product_type='DpdLE3clAmicoAux')
  >>> print("Found", len(le3_product_list), "results")
  Found 2 results
  >>> print(le3_product_list)
  basic_download_data_oid   product_type                      product_id                    observation_id_list tile_index_list patch_id_list filter_name
  ----------------------- ---------------- ------------------------------------------------ ------------------- --------------- ------------- -----------
                    47257 DpdLE3clAmicoAux PPO_REGREPROC1_R2_CLTEST_R0_CLDET_R3-amico_aux-0                  {}              {}            {}
                    47258 DpdLE3clAmicoAux PPO_REGREPROC1_R2_CLTEST_R0_CLDET_R7-amico_aux-0                  {}              {}            {}









2. Authenticated access
-----------------------

Authenticated users are able to access to TAP+ capabilities (shared tables, persistent jobs, etc.) In order to
authenticate a user, ``login`` method must be called. After a successful authentication, the user will be authenticated
until the ``logout`` method is called.

All previous methods (``query_object``, ``cone_search``, ``load_table``, ``load_tables``, ``launch_job``) explained for
non authenticated users are applicable for authenticated ones.

The main differences are:

* Asynchronous results are kept at the server side forever (until the user decides to remove one of them).
* Users can access to shared tables.


2.1. Login/Logout
^^^^^^^^^^^^^^^^^

There are several ways to log in to the Euclid archive.

**Login through graphic interface**

*Note: The Python Tkinter module is required to use the login_gui method.*

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> from astroquery.esa.euclid import Euclid
  >>> Euclid.login_gui()


**Login through command line**

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login()
  >>> User: user
  >>> Password: pwd (not visible)

or

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login(user='userName', password='userPassword')


It is possible to use a file where the credentials are stored:

*The file must contain user and password in two different lines.*

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login(credentials_file='my_credentials_file')

To perform a logout:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.logout()



2.2. User space management: table upload
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is now possible to store a table in the private user space. The table to be uploaded can
be in a VOTable_ located at a given URL, a table stored in a local file in the user machine,
a pre-computed Astropy table file or a job executed in the Euclid archive.


Each user has a database schema described as: 'user_<user_login_name>'. For instance, if a
login name is 'joe', the database schema is 'user_joe'. Your uploaded table can be
referenced as 'user_joe.table_name'

2.2.1. Uploading table from URL
"""""""""""""""""""""""""""""""""""

An already generated VOTable, accessible through a URL, can be uploaded to Euclid archive.

The following example launches a query to Vizier TAP ('url' parameter). The result is a
VOTable that can be uploaded to the user's private area.

Your schema name will be automatically added to the provided table name:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login()
  >>> # Provide a URL pointing to valid VOTable resource
  >>> url = ("https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync/?"
  ...        "REQUEST=doQuery&lang=ADQL&FORMAT=votable&"
  ...        "QUERY=select+*+from+TAP_SCHEMA.columns+where+table_name='II/336/apass9'")
  >>> job = Euclid.upload_table(upload_resource=url, table_name="table_test_from_url",
  ... table_description="Some description")
  Job '1539932326689O' created to upload table 'table_test_from_url'.

Now, you can query your table as follows (a full qualified table name must be provided,
i.e.: *user_<your_login_name>.<table_name>*. Note that if the <table_name> contains capital letters, it must be
surrounded by quotation marks, i.e.: *user_<your_login_name>."<table_name>"*):

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> full_qualified_table_name = 'user_<your_login_name>.table_test_from_url'
  >>> query = 'select * from ' + full_qualified_table_name
  >>> job = Euclid.launch_job(query=query)
  >>> results = job.get_results()


2.2.2. Uploading table from file
"""""""""""""""""""""""""""""""""""

A file containing a table can be uploaded to the user private area. Only a file associated to any of the formats described in
https://docs.astropy.org/en/stable/io/unified.html#built-in-table-readers-writers, and automatically identified by its suffix
or content can be used. Note that for a multi-extension fits file with multiple tables, the first table found will be used.
For any other format, the file can be transformed into an astropy Table (https://docs.astropy.org/en/stable/io/unified.html#getting-started-with-table-i-o)
and passed to the method.

The parameter 'format' must be provided when the input file is not a votable file.

Your schema name will be automatically added to the provided table name.

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login()
  >>> job = Euclid.upload_table(upload_resource="1535553556177O-result.vot", table_name="table_test_from_file", format="votable")
  Sending file: 1535553556177O-result.vot
  Uploaded table 'table_test_from_file'.

Now, you can query your table as follows (a full qualified table name must be provided,
i.e.: *user_<your_login_name>.<table_name>*. Note that if the <table_name> contains capital letters, it must be
surrounded by quotation marks, i.e.: *user_<your_login_name>."<table_name>"*):

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> full_qualified_table_name = 'user_<your_login_name>.table_test_from_file'
  >>> query = 'select * from ' + full_qualified_table_name
  >>> job = Euclid.launch_job(query=query)
  >>> results = job.get_results()


2.2.3. Uploading table from an astropy Table
"""""""""""""""""""""""""""""""""""

A votable can be uploaded to the server in order to be used in a query. Your schema name will be automatically added to the provided table name.

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> from astropy.table import Table
  >>> a=[1,2,3]
  >>> b=['a','b','c']
  >>> table = Table([a,b], names=['col1','col2'], meta={'meta':'first table'})
  >>> # Upload
  >>> Euclid.login()
  >>> Euclid.upload_table(upload_resource=table, table_name='table_test_from_astropy')


Now, you can query your table as follows (a full qualified table name must be provided,
i.e.: *user_<your_login_name>.<table_name>*. Note that if the <table_name> contains capital letters, it must be
surrounded by quotation marks, i.e.: *user_<your_login_name>."<table_name>"*):

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> full_qualified_table_name = 'user_<your_login_name>.table_test_from_astropy'
  >>> query = 'select * from ' + full_qualified_table_name
  >>> job = Euclid.launch_job(query=query)
  >>> results = job.get_results()


2.2.4. Uploading table from job
"""""""""""""""""""""""""""""""""""

The results generated by an *asynchronous* job (from a query executed in the Euclid archive) can be
ingested in a table in the user's private area.

The following example generates a job in the Euclid archive and then, the results are ingested in a
table named: user_<your_login_name>.'t'<job_id>:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login()
  >>> job_1 = Euclid.launch_job_async("select top 10 * from catalogue.mer_catalogue")
  >>> Euclid.upload_table_from_job(job=job_1)
  Created table 't1539932994481O' from job: '1539932994481O'.

Now, you can query your table as follows (a full qualified table name must be provided,
i.e.: *user_<your_login_name>."t<job_id>"*. Note that the previous table name must be
surrounded by quotation marks since it contains capital letters.):

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> full_qualified_table_name = 'user_<your_login_name>."t1710251325268O"'
  >>> query = 'select * from ' + full_qualified_table_name
  >>> job = Euclid.launch_job(query=query)
  >>> results = job.get_results()


2.2.5 Deleting table
"""""""""""""""""""""""""""""""""""

A table from the user's private area can be deleted as follows:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login_gui()
  >>> job = Euclid.delete_user_table(table_name="table_test_from_file")
  Table 'table_test_from_file' deleted.

2.2.6 Updating table metadata
"""""""""""""""""""""""""""""""""""

It can be useful for the user to modify the metadata of a given table. For example, a user
might want to change the description (UCD) of a column, or the flags that give extra information
about a certain column. This is possible using:

.. Skipping authentication requiring examples
.. doctest-skip::
  >>> Euclid.login_gui()
  >>> Euclid.update_user_table(table_name, list_of_changes)

where the list of changes is a list of 3 items:

["column name to be changed", "metadata parameter to be changed", "new value"]

The metadata parameter to be changed can be 'utype', 'ucd', 'flags' or 'indexed':

* values for 'utype' and 'ucd' are free text. See VOTable_ specification (sections UType and UCD), UCD_ specification and UTypes_ usage.

* value for 'flags' can be 'Ra', 'Dec', 'Mag', 'Flux' and 'PK'.

* value for 'indexed' is a boolean indicating whether the column is indexed or not.

.. _UCD: https://www.ivoa.net/documents/latest/UCD.html
.. _UTypes: https://www.ivoa.net/documents/Notes/UTypesUsage/index.html


It is possible to apply multiple changes at once.
This is done by putting each of the changes in a list. See example below.

In this case, we have a table (user_joe.table), with several columns: 'recno', 'nobs', 'raj2000' and 'dej2000'.

We want to set:

* 'ucd' of 'recno' column to 'ucd sample'
* 'utype' of 'nobs' column to 'utype sample'
* 'flags' of 'raj2000' column to 'Ra'
* 'flags' of 'dej2000' column to 'Dec'

We can type the following:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login_gui()
  >>> Euclid.update_user_table(table_name="user_joe.table",
  ...                        list_of_changes=[["recno", "ucd", "ucd sample"],
  ...                                         ["nobs","utype","utype sample"],
  ...                                         ["raj2000","flags","Ra"],
  ...                                         ["dej2000","flags","Dec"]])
  Retrieving table 'user_joe.table'
  Parsing table 'user_joe.table'...
  Done.
  Table 'user_joe.table' updated.


2.3. Tables sharing
^^^^^^^^^^^^^^^^^^^

It is possible to share tables with other users. You have to create a group, populate that
group with users, and share your table to that group. Then, any user belonging to that group
will be able to access your shared table in a query.

2.3.1. Creating a group
"""""""""""""""""""""""""""""""""""

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login()
  >>> Euclid.share_group_create(group_name="my_group", description="description")

2.3.2. Removing a group
"""""""""""""""""""""""""""""""""""

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.share_group_delete(group_name="my_group")

2.3.3. Listing groups
"""""""""""""""""""""""""""""""""""

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> groups = Euclid.load_groups()
  >>> for group in groups:
  ...     print(group.title)

2.3.4. Adding users to a group
"""""""""""""""""""""""""""""""""""

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.share_group_add_user(group_name="my_group",user_id="<user_login_name")

2.3.5. Removing users from a group
"""""""""""""""""""""""""""""""""""

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.share_group_delete_user(group_name="my_group",user_id="<user_login_name>")

2.3.6. Sharing a table to a group
"""""""""""""""""""""""""""""""""""

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.share_table(group_name="my_group",
  ...                  table_name="user_<user_login_name>.my_table",
  ...                  description="description")

2.3.7. Stop sharing a table
"""""""""""""""""""""""""""""""""""

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.share_table_stop(table_name="user_<user_login_name>.my_table", group_name="my_group")


2.3.8 Listing shared tables
"""""""""""""""""""""""""""""""""""

In the Euclid archive user tables can be shared among user groups.

To obtain a list of the tables shared to a user type the following:

  >>> tables = Euclid.load_tables(only_names=True, include_shared_tables=True)
  >>> for table in tables:
  ...   print(table.get_qualified_name())
    catalogue.mer_catalogue
    catalogue.mer_cutouts
    catalogue.mer_morphology
    catalogue.phz_classification
    catalogue.phz_galaxy_sed
    ...

.. _uploading_table_to_user_space:




2.5. Cross match
^^^^^^^^^^^^^^^^

It is possible to run a geometric cross-match between the RA/Dec coordinates of two tables using the crossmatch function
provided by the archive. The returned table includes the identifiers from both tables and the angular separation, in
degrees, between the RA/Dec coordinates of each source in the first table and its corresponding source in the second
table.

The cross-match can be executed in one single step by the following method

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login()
  >>> full_qualified_table_name = 'user_<your_login_name>.my_sources'
  >>> job = Euclid.cross_match_basic(table_a_full_qualified_name=full_qualified_table_name, table_a_column_ra='raj2000',
                     table_a_column_dec='dej2000', table_b_full_qualified_name='catalogue.mer_catalogue',
                     table_b_column_ra='right_ascension', table_b_column_dec='declination, radius=1.0, background=True)
  >>> result = job.get_results()

This method updates the user table metadata to flag the positional RA/Dec columns and launches the positional
cross-match as an asynchronous query. The returned job provides direct access to the output of the cross-match
information: for each matched source, all the columns from the input tables plus the angular distance (degrees).
Therefore, the size of the output can be quite large.

By default, this method targets the main catalogue associated to each environment (PDR, OTF, REG and IDR) using a cone
search radius of 1.0 arcseconds. Therefore, the above example can also be simplified as follows

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> Euclid.login()
  >>> full_qualified_table_name = 'user_<your_login_name>.my_sources'
  >>> job = Euclid.cross_match_basic(table_a_full_qualified_name=full_qualified_table_name, table_a_column_ra='raj2000',
                                   table_a_column_dec='dej2000')
  >>> result = job.get_results()


3. DataLink service
-----------------------------------

DataLink_ is a data access protocol compliant with the IVOA_ architecture that provides a linking mechanism between
datasets offered by different services. In practice, it can be seen and used as a web service providing the list of additional
data products available for each object outside the main catalogue(s). For more information about the products served via
DataLink in the Euclid ESA Archive we recommend reading the Archive DataLink tutorials available at https://eas.esac.esa.int/sas/.

The DataLink products are restricted to authenticated users via the `~astroquery.utils.tap.TapPlus.load_data` method.
From SAS the Datalink service can be used to access and download 1D Spectra data.

To find out the resources associated with a given source:

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> ids=["2707008224650763513"]
  >>> datalink = Euclid.get_datalinks(ids=ids)
  >>> print(datalink)
             ID            linking_parameter                                          access_url                                         service_def error_message semantics     description     content_type content_length
                                                                                                                                                                                                                   byte
  ------------------------ ----------------- ------------------------------------------------------------------------------------------- ----------- ------------- --------- ------------------- ------------ --------------
  sedm 2707008224650763513         SOURCE_ID https://eas.esac.esa.int/sas-dd/data?ID=sedm+2707008224650763513&RETRIEVAL_TYPE=SPECTRA_RGS                               #this  Spectra Red Source                          --
  sedm 2707008224650763513         SOURCE_ID https://eas.esac.esa.int/sas-dd/data?ID=sedm+2707008224650763513&RETRIEVAL_TYPE=SPECTRA_BGS                               #this Spectra Blue Source                          --




The query below retrieves a random sample of Euclid sources having spectra.

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> from astroquery.esa.euclid import Euclid
  >>> query = f"SELECT TOP 2000 * FROM catalogue.spectra_source"
  >>> job = Euclid.launch_job_async(query)
  >>> results = job.get_results()
  >>> print(f'Table size (rows): {len(results)}')
  Table size (rows): 2000
  >>> print(results)
  combined_spectra_fk combined_spectra_product_fk            datalabs_path                 dec_obj      dith_num                           file_name                           ... hdu_index      ra_obj           source_id      spectra_source_oid to_be_published
  ------------------- --------------------------- ----------------------------------- ----------------- -------- ------------------------------------------------------------- ... --------- ---------------- ------------------- ------------------ ---------------
                  161                        6170 /data/euclid_q1/Q1_R1/SIR/102159190  66.2272115578693        3 EUC_SIR_W-COMBSPEC_102159190_2024-11-05T16:27:45.906227Z.fits ...       355 267.261146414289 2672611464662272115              66161               1
                  161                        6170 /data/euclid_q1/Q1_R1/SIR/102159190   66.230432248046        4 EUC_SIR_W-COMBSPEC_102159190_2024-11-05T16:27:45.906227Z.fits ...       557 267.319331443563 2673193314662304322              66179               1
                  161                        6170 /data/euclid_q1/Q1_R1/SIR/102159190  66.2259968885041        4 EUC_SIR_W-COMBSPEC_102159190_2024-11-05T16:27:45.906227Z.fits ...       679  267.39974379438 2673997437662259968              66185               1
                  ...                         ...                                 ...               ...      ...                                                           ... ...       ...              ...                 ...                ...             ...
  Length = 2000 rows
  >>> print("source ids:")
  >>> print(results['source_id'])
  <Column name='source_id' dtype='int64' length=200>
  2672611464662272115
  2673193314662304322
  2673997437662259968
                  ...


The following example shows how to retrieve the DataLink products (1D Spectra) associated with the previous sources (IDs).

.. Skipping authentication requiring examples
.. doctest-skip::

  >>> files = Euclid.get_spectrum(retrieval_type='SPECTRA_BGS', source_id='2675005060662306333')
  >>> from astropy.io import fits
  >>> print(fits.info(files[0]))  # doctest: +IGNORE_OUTPUT
  Filename: /home/astroquery/temp_20250225_204959/2675005060662306333.fits
  No.    Name      Ver    Type      Cards   Dimensions   Format
    0  PRIMARY       1 PrimaryHDU       4   ()
    1  2675005060662306333    1 BinTableHDU     35   531R x 6C   [1E, 1E, 1J, 1E, 1E, 1I]
  None


A fits file is made if no file name is provided.


.. _appendix:

========
Appendix
========

The following table summarises the available values of the parameters of the method get_scientific_product_list.

.. csv-table:: Valid values for the parameters of the method get_scientific_product_list
    :file: table_values.csv
    :header-rows: 1
    :widths: auto

=============
Reference/API
=============

.. automodapi:: astroquery.esa.euclid
    :no-inheritance-diagram:
