pro RIPS_Mercury_Perkins_v01, PART=part, NIGHT=night

; *********************************************************************************************************************
; *********************************************************************************************************************
; Quick routine to extract Na emission from spectral scans over Mercury's disk with Perkins/RIPS in order to create a 2-D
; image of Na above and near Mercury's disk.
; 
; The program is split into 4 "parts":
;   0 = break kinetic series into manageable datacube
;       input:   kinetic series FITS file(s)
;       output:  RIPS imaging (x,y,t) and spectral (wl,y,t) datacubes 
;   1 = find Mercury centroids by cross-corellation with previous image
;       input:   imaging datacube from Part 0
;       output:  co-aligned imaging datacube (x,y,t)
;   2 = isolate sodium emission in every frame
;       input:   spectral datacube from Part 0
;       output:  the calibrated spectral datacube (wl,y,t) for Na exosphere
;   3 = extract Na signal and place into 1D spectra
;       input:   spectral datacube from Part 2 (wl,y,t)
;       output:  Na brightness (wl,y) and linewidth (wl,y)
;   4 = put everything together and build exosphere image
;       input:   Na brightness and linewidth from Part 3 and imaging cube from Part 0
;       output:  none yet, a combined image called "test" for now...
;  99 = do all of 0-4
;
; Input variables
;   part           - the "part" of the program to execute
;   night          - the night of AEOS data to use; either '20' or '25' for June 20 or 25 (Default='25'); 
;                    contains numerous sub-variables: mercury_dir, dark_dir, flat_dir, sky_dir, Mercury_file, 
;                    sky_file, dark_file, flat_file, Na_D_rngs, ims, sps
;   thresh         - hotpixel removal threshold (sigma above background; Default=12)
;   width          - median smoothing width to get local background reference (for pixels above thresh; Default=5)
;
; 
; Created:  18 Oct 2018 - copied from RIPS_Mercury_AEOS_v03.pro to start
; Edits:    
; *********************************************************************************************************************
; *********************************************************************************************************************
; =====================================================================================================================
; Define variables
; =====================================================================================================================
SetDefaultValue, part, 2
SetDefaultValue, night, '15'                                              ; set to '20' or '25' for night of June 20 or June 25
SetDefaultValue, thresh, 12                                               ; hotpixel removal threshold (sigma above background)
SetDefaultValue, width, 5                                                 ; Median smoothing width to get local background reference (for pixels above thresh) 
SetDefaultValue, minfac, 0.8                                              ; maximum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, maxfac, 1.2                                              ; minimum factor to apply to reference spectrum for Na emission extraction
SetDefaultValue, dfac, 0.01                                               ; factor increment for minfac->maxfac
SetDefaultValue, ct, 20;27                                                   ; default color table
SetDefaultValue, max_spec, 2500.;2950.                                          ; partial hack --> the max value of the spectral cube after coadding (used for "prettier" movies)
SetDefaultValue, max_img, 1.6e6;2.7e6                                           ; partial hack --> the max value of the image cube after coadding (used for "prettier" movies)
;SetDefaultValue, img_extraction, [275,425,75,225]                        ; coordinates for extracting a square Mercury image, [x1,x2,y1,y2]
;SetDefaultValue, img_extraction, [165,315,85,235]                        ; coordinates for extracting a square Mercury image, [x1,x2,y1,y2]
SetDefaultValue, img_extraction, [145,295,85,235]                        ; coordinates for extracting a square Mercury image, [x1,x2,y1,y2]
;SetDefaultValue, img_extraction, [65,164,51,150]                        ; coordinates for extracting a square Mercury image, [x1,x2,y1,y2] -- attempt at smaller cutouts to speed things up
;SetDefaultValue, spec_extraction, [156,316,120,280]                       ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
;SetDefaultValue, spec_extraction, [156,316,80,240]                       ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
;SetDefaultValue, spec_extraction, [165,315,85,235]                       ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
;SetDefaultValue, spec_extraction, [145,295,85,235]                       ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
SetDefaultValue, spec_extraction, [145,295,77,227]                       ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2]
;SetDefaultValue, spec_extraction, [65,164,51,150]                       ; coordinates for extracting a square Mercury spectral image, [x1,x2,y1,y2] -- attempt at smaller cutouts to speed things up
SetDefaultValue, movie_scale, 4                                           ; scale the x-y dimensions of the extracted "Na image" and "Mercury disk" arrays for movies by this value
SetDefaultValue, do_realign, 0                                            ; 0=run part 1 as normal, 1=use previous cross-correlation as master template
SetDefaultValue, do_rotate, 0                                             ; 0=no rotation of image, 1=rotate Na image (inverted, =5 -- see below)
SetDefaultValue, do_median, 0                                             ; 0=use mean when constructing Na Mercury image, 1=use median instead
SetDefaultValue, do_smooth, 1                                             ; set to 1 to smooth reference and frame images before coaligning
SetDefaultValue, stddev_cutoff, 230.                                      ; only use Mercury frames with stddev > stddev_cutoff, these are defined as "good" 

mercury_dir   = 'Y:\obs_18\Perkins_RIPS_March\15\renamed\'                ; directory with Mercury FITS files
dark_dir      = 'Y:\obs_18\Perkins_RIPS_March\14\Carl_Keep\'              ; directory with dark file
flat_dir      = 'Y:\obs_18\Perkins_RIPS_March\15\'                        ; directory with flat file
sky_dir       = 'Y:\obs_18\Perkins_RIPS_March\15\'                        ; directory with sky file
;Mercury_file  = 'RIPS1_Thu Mar 15 2018_01.49.16_791.fits'                 ; Mercury kinetic series frames (see "March observing logs.xlsx" for more info)
Mercury_file  = 'Mercury_Na_' + strfix(indgen(16)+789) + '.fits'   ; this line for "RIPS1...789.fits" to "...804.fits"
Mercury_file  = 'Mercury_Na_' + strfix(indgen(5)+789) + '.fits'
Mercury_file  = 'Mercury_Na_' + strfix(['791','792','793','794','798','800','803']) + '.fits'
;Mercury_file  = 'Mercury_Na_' + strfix(['791','792']) + '.fits'
sky_file      = 'RIPS1_Thu Mar 15 2018_01.12.12_785.fits'                 ; sky frame
dark_file     = 'Dark.fits'                                               ; dark frame
flat_file     = 'Flat_NaSpectra_624slitwidth_NaND1imaging1.fits'          ; flat frame
;Na_D_rngs     = [266,297,528,538]                                         ; ranges of Na D Mercury emission (x pixels in spec domain for D1 and D2)
Na_D_rngs     = [235,255,475,495]
;Na_D_rngs     = [240,250,240,250]
frame_ref     = 182                                                       ; this number identifies a particular frame of imaging_cube as the reference frame (used for aligning all other frames spatially)
offset_maxes  = [80,40]                                                   ; if the absolute value of the [x,y] offset is larger than these values then that frame is marked "bad"
ims           = [283,759,39,374]                                          ; variables for image analysis: x1,x2,y1,y2 (formerly "imaging_statsec")
;ims           = [383,659,89,324]                                          ; variables for image analysis: x1,x2,y1,y2 (formerly "imaging_statsec") -- attempt at smaller cutouts to speed things up
sps           = [99,899,502,941]                                          ; variables for spectral analysis: x1,x2,y1,y2 (formerly "spectra_statsec")
;sps           = [199,799,552,891]                                          ; variables for spectral analysis: x1,x2,y1,y2 (formerly "spectra_statsec") -- attempt at smaller cutouts to speed things up
iref_factor   = 1.                                                        ; factor to scale reference spectrum by
count_factor  = 1.e4                                                      ; this is the scaling factor to apply to the white light Mercury image (i.e., the slit covers some regions of the disk, so we scale based on total slit coverage at each location)
integration   = 0.05                                                      ; HACK (in that we should just look up the integration time from the header... but I'm lazy)
n_frames_per_file = 100                                                   ; # of frames for each "Mercury_file" from above, typically 100

outdir        =  Mercury_dir + 'Processed\'                               ; output directory
nfiles        =  n_elements(Mercury_file)                                 ; # of different Mercury kinetic frames

!p.multi = 0.
cgloadct, ct;, /reverse

; =====================================================================================================================
; Find display size, set windows accordingly
; =====================================================================================================================
DEVICE, GET_SCREEN_SIZE = ss ;& PRINT, ss                                 ; find screen size
winpos_x      = [ss(0)/2, ss(0)/2, ss(0)/2+ss(0)/4, ss(0)/2+ss(0)/4, ss(0)] ; default x window locations
;winpos_x      = [ss(0)/3, ss(0)/3, ss(0)/3+ss(0)/4, ss(0)/3+ss(0)/4, ss(0)] ; default x window locations
winpos_y      = [0, ss(1)/2, 0, ss(1)/2, 0]                               ; default y window locations

; =====================================================================================================================
; Part 0 : break kinetic series into manageable datacube
; =====================================================================================================================
if part eq 0 or part eq 99 then begin 
    icube = MRDFITS(Mercury_dir + Mercury_file(0), 0, header, /unsigned, /silent )
    s = size(icube, /dimensions) 
    cube = fltarr(s(0),s(1),s(2)*nfiles)
    for ifile = 0, nfiles-1 do begin
      icube = MRDFITS(Mercury_dir + Mercury_file(ifile), 0, header, /unsigned, /silent )
      cube(*,*,s(2)*ifile:s(2)*ifile+n_frames_per_file-1) = icube
    endfor
    
    imaging_cube = reform(cube[ims(0):ims(1),ims(2):ims(3), *])
    spectra_cube = reform(cube[sps(0):sps(1),sps(2):sps(3), *])
    MWRFITS, imaging_cube, outdir + 'imaging_cube.fits', header, /create, /silent
    MWRFITS, spectra_cube, outdir + 'spectra_cube.fits', header, /create, /silent
    Sky_cube = MRDFITS(sky_dir + sky_file, 0, header, /unsigned, /silent )
    imaging_sky_cube = reform(sky_cube[ims(0):ims(1),ims(2):ims(3), *])
    spectra_sky_cube = reform(sky_cube[sps(0):sps(1),sps(2):sps(3), *])
    MWRFITS, imaging_sky_cube, outdir + 'imaging_sky_cube.fits', header, /create, /silent
    MWRFITS, spectra_sky_cube, outdir + 'spectra_sky_cube.fits', header, /create, /silent
    beep
;    stop
endif

; =====================================================================================================================
; Part 1 : find Mercury centroids by cross-corellation with previous image
; =====================================================================================================================
if part eq 1 or part eq 99 then begin ;Find the centroids by cross-correlation with the last image
    if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
    if part ne 99 then spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
    frame = reform(imaging_cube[*,*,0])
    s = size(imaging_cube, /dimensions)
    ss = size(spectra_cube, /dimensions) 
    window, 0, xs = s[0], ys = s[1], xpos=winpos_x[0],         ypos=winpos_y[0]
    window, 1, xs = s[0], ys = s[1], xpos=winpos_x[0],         ypos=winpos_y[0]+s[1]+40
    window, 2, xs = s[0], ys = s[1], xpos=winpos_x[0]+s[0]+20, ypos=0
    window, 3, xs = s[0], ys = s[1], xpos=winpos_x[0]+s[0]+20, ypos=s[1]+40

    iDark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )
    dark = idark[ims(0):ims(1),ims(2):ims(3)]
    specdark = idark[sps(0):sps(1),sps(2):sps(3)]
    dark = sigma_filter( dark, 5, N_sigma=5)
    specdark = sigma_filter( specdark, 5, n_sigma=5 )
    iflat = MRDFITS(flat_dir + flat_file, 0, header, /fscale ) 
    ; flat = (flat[ims(0):ims(1),ims(2):ims(3)] - dark) / mean(flat[ims(0):ims(1),ims(2):ims(3)])
    ; Above line was from RIPS_Mercury_AEOS..., removed for now in favor of the following 2 lines
    flat = ( iflat[ims(0):ims(1),ims(2):ims(3)] - dark ) / max( iflat[ims(0):ims(1),ims(2):ims(3)] - dark )
    specflat = ( iflat[sps(0):sps(1),sps(2):sps(3)] - specdark ) / max( iflat[sps(0):sps(1),sps(2):sps(3)] - specdark )
    flat = flat + (1. - median(flat))            ; this normalizes such that the median of the flat is 1
    specflat = specflat + (1. - median(specflat))
    flat[WHERE(flat lt .01, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
    flat[WHERE(flat gt 2.0, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
    specflat[WHERE(specflat lt .01, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
    specflat[WHERE(specflat gt 2.0, /NULL)] = !values.f_Nan ;reject unusual counts for centroid
;    flat = smart_shift(flat, .8, .5,  /interp) ;hack
    gooddata = where(Finite(flat), ngooddata, comp=baddata, ncomp=nbaddata)
    if nbaddata gt 0 then flat[baddata] = interpol(flat[gooddata], gooddata, baddata)
    gooddata_spec = where(Finite(specflat), ngooddata_spec, comp=baddata_spec, ncomp=nbaddata_spec)
    if nbaddata_spec gt 0 then specflat[baddata_spec] = interpol(specflat[gooddata_spec], gooddata_spec, baddata_spec)
    if do_realign then begin
      restore, outdir+'shift_array.sav'          ; contains shift_array, std_devs, and aligned_imaging_cube
      shift_array_old = shift_array              ; keep old shift_array for comparison
      reference = total(aligned_imaging_cube,3)  ; create "master Mercury" template for cross-correlation
    endif else begin
;      ireference = (reform(imaging_cube[*,*,s[2]-1]) - dark) / flat
      ireference = (reform(imaging_cube[*,*,frame_ref]) - dark) / flat
      acre, ireference, reference, thresh, width
    endelse
    
    slit_width = 5                                              ; hack -- width of slit
    xslit = where( total(flat,2) eq min(total(flat,2)) )        ; the "middle" of the slit

    wset, 1
    cgimage, reference, /axes, title='Reference image for alignment'
    shift_array = intarr(s[2],2)
    std_devs = fltarr(s[2])
    for i = 0, s[2]-1 do begin
        iframe = (reform(imaging_cube[*,*,i]) - dark) / flat
        iframe[xslit-slit_width:xslit+slit_width,*] = !values.F_Nan
        fill_missing, iframe, !values.F_Nan, 1                      ; interpolate over slit
        ispecframe = (reform(spectra_cube[*,*,i]) - specdark) / specflat
        acre, iframe, frame, thresh, width
        acre, ispecframe, specframe, thresh, width
;        CORREL_OPTIMIZE, reference, frame, xoffset_optimum, yoffset_optimum, NUMPIX = 150
        if do_smooth then begin
          frame_bandpass = bandpass_filter(smooth(frame,8,/edge_truncate), 0., 0.15, /butterworth)
          ref_bandpass = bandpass_filter(smooth(reference,8,/edge_truncate), 0., 0.15, /butterworth)
;          kernel = gaussian_function([3,3], width=55, maximum=255)
;          ref_gauss = convol(ref_bandpass,kernel,invalid=255,missing=0,/edge_zero,/normalize)
;          frame_gauss = convol(frame_bandpass,kernel,invalid=255,missing=0,/edge_zero,/normalize)
        endif else begin
          frame_bandpass = bandpass_filter(frame, 0., 0.15, /butterworth)
          ref_bandpass = bandpass_filter(reference, 0., 0.15, /butterworth)
;          stop  ; stop for now because I'm testing the smoothing (LM)
        endelse
; next 4 lines TEMP HACK!
;ref_gauss = reference
;frame_gauss = frame
;ref_bandpass = ref_gauss
;frame_bandpass = frame_gauss
;        CORREL_OPTIMIZE, reference^2., frame^2., xoffset_optimum, yoffset_optimum, NUMPIX=150   ; let's try squaring them...
        CORREL_OPTIMIZE, ref_bandpass, frame_bandpass, xoffset_optimum, yoffset_optimum, NUMPIX=150  ; let's try squaring them...
;        CORREL_OPTIMIZE, ref_bandpass^2., frame_bandpass^2., xoffset_optimum, yoffset_optimum, NUMPIX=150  ; let's try squaring them...
;        CORREL_OPTIMIZE, ref_gauss, frame_gauss, xoffset_optimum, yoffset_optimum, NUMPIX=150  ; let's try squaring them...
;        if abs(xoffset_optimum) gt offset_maxes[0] or abs(yoffset_optimum) gt offset_maxes[1] then shift_array[i,*] = shift_array[i-1,*]
        shift_array[i,*] = [xoffset_optimum, yoffset_optimum]
        if abs(xoffset_optimum) gt offset_maxes[0] or abs(yoffset_optimum) gt offset_maxes[1] then begin;shift_array[i,*] = [-666,-666]  ; ie this is a "bad" frame
          xoffset_optimum = (where(total(ref_bandpass^2.,2) eq max(total(ref_bandpass^2.,2))))[0] $
                          - (where(total(frame_bandpass^2.,2) eq max(total(frame_bandpass^2.,2))))[0]
          yoffset_optimum = (where(total(ref_bandpass^2.,1) eq max(total(ref_bandpass^2.,1))))[0] $
                          - (where(total(frame_bandpass^2.,1) eq max(total(frame_bandpass^2.,1))))[0]
          shift_array[i,*] = [xoffset_optimum, yoffset_optimum]
          print, 'CORREL failed for frame ' + strfix(i) + ' -- new offsets (x,y): ' + strfix(xoffset_optimum) + ', ' + strfix(yoffset_optimum)
        endif else print, i, shift_array(i,0), shift_array(i,1)
;        print, i, shift_array(i,0), shift_array(i,1)
        
        wset, 0
;        cgimage, smooth(ref_bandpass(180:250,110:200),8), /axes, title='initial comparison'
;        cgimage, ref_bandpass(180:250,110:200)/max(ref_bandpass), /axes, title='initial bandpass comparison'
        cgcontour, ref_bandpass(180:250,110:200)/max(ref_bandpass), title='initial bandpass comparison', levels=[0.2,0.4,0.6,0.8], /fill
        cgcontour, ref_bandpass(180:250,110:200)/max(ref_bandpass), title='initial bandpass comparison', levels=[0.2,0.4,0.6,0.8], /onimage
        cgcontour, frame_bandpass(180:250,110:200)/max(frame_bandpass), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='red'
;        cgcontour, smooth(frame_bandpass(180:250,110:200),8,/edge_truncate), /onimage
        wset, 2
;        cgimage, smooth(ref_bandpass(180:250,110:200),8), /axes, title='aligned'
;        cgimage, ref_bandpass(180:250,110:200)/max(ref_bandpass), /axes, title='aligned bandpass images', levels=[0.2,0.4,0.6,0.8]
        cgcontour, ref_bandpass(180:250,110:200)/max(ref_bandpass), title='aligned bandpass images', levels=[0.2,0.4,0.6,0.8], /fill
        cgcontour, ref_bandpass(180:250,110:200)/max(ref_bandpass), title='aligned bandpass images', levels=[0.2,0.4,0.6,0.8], /onimage
        f = shift(frame_bandpass, [shift_array[i,*]])
;        cgcontour, smooth(f(180:250,110:200),8,/edge_truncate), /onimage
        cgcontour, f(180:250,110:200)/max(f), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='red'
        wset, 1
        cgimage, reference(180:250,110:200), /axes, title='raw reference'
        wset, 3
        f = shift(frame, [shift_array[i,*]])
        cgimage, f(180:250,110:200), /axes, title='aligned raw frame'
        cgcontour, ref_bandpass(180:250,110:200)/max(ref_bandpass), levels=[0.2,0.4,0.6,0.8], /onimage
;        stop
    endfor

    good_frames = where(shift_array(*,0) ne -666)
    aligned_imaging_cube = fltarr(s[0],s[1],n_elements(good_frames))      ; imaging cube after everything is aligned
    aligned_spectra_cube = fltarr(ss[0],ss[1],n_elements(good_frames))    ; spectral cube after realignment (note that +y in the spatial dimension is -y in the spectral!)
    aligned_imaging_cube2 = fltarr(s[0],s[1],n_elements(good_frames))     ; another version, this time with the slit pixels zeroed out

    for i = 0, n_elements(good_frames)-1 do begin
        frame = (reform(imaging_cube[*,*,i]) - dark) / flat ;make this into a movie!
        specframe = (reform(spectra_cube[*,*,i])); - specdark) / specflat
        aligned_imaging_cube[*,*,i] = shift(frame, [shift_array[i,*]])
        aligned_spectra_cube[*,*,i] = shift(specframe, [0,-shift_array[i,1]])  ; note that we are taking the negative of the y shift from the image!
;        wset, 0
;        cgimage, total(aligned_imaging_cube,3), /axes, title='Total of aligned frames: 0->' + strfix(i)
        img_frame = reform(frame(img_extraction(0)-50:img_extraction(1)+50,img_extraction(2)-50:img_extraction(3)+50))   ; extract just the image part of the frame
        aligned_imaging_cube2[*,*,i] = shift(frame, [shift_array[i,*]])
        std_devs(i) = stddev(img_frame)                                   ; seems like anything < 1000 is not that good
;        print, i, n_elements(frame) - n_elements( where( frame lt (max(frame)*.8)/2d ) ), shift_array(i,0), shift_array(i,1)
;        wset, 2
;        cgimage, hist_equal(frame), /axes, title='Raw frame: ' + strfix(i)
;        pmax = where(frame eq max(frame))
;        xmax = ( where_2D(frame, pmax) )[0]
;        ymax = ( where_2D(frame, pmax) )[1]
;        cgoplot, [1,1]*xmax, !Y.Crange, col='yellow', thick=2
;        cgoplot, !X.Crange, [1,1]*ymax, col='yellow', thick=2
 ;       cgtext, 5, ymax+5, strfix(xmax), charsize=1.5
 ;       cgtext, xmax+5, 5, strfix(ymax), charsize=1.5
 ;       wset, 3
;        aframe = shift(frame, [shift_array[i,*]])
;        cgimage, hist_equal(aframe), /axes, title='Aligned frame: ' + strfix(i)
;        pmax = where(aframe eq max(aframe))
;        xmax = ( where_2D(aframe, pmax) )[0]
 ;       ymax = ( where_2D(aframe, pmax) )[1]
 ;       cgoplot, [1,1]*xmax, !Y.Crange, col='yellow', thick=2
 ;       cgoplot, !X.Crange, [1,1]*ymax, col='yellow', thick=2
 ;       cgtext, 5, ymax+5, strfix(xmax), charsize=1.5
 ;       cgtext, xmax+5, 5, strfix(ymax), charsize=1.5
 ;       wait, 0.15
;        stop
    endfor    
    save, shift_array, aligned_imaging_cube, aligned_spectra_cube, std_devs, filename = outdir + 'shift_array.sav'
    beep
;    stop
;    x = reform(aligned_imaging_cube(img_extraction(0)-50:img_extraction(1)+50,img_extraction(2)-50:img_extraction(3)+50,i))
endif

; =====================================================================================================================
; Part 2 : Isolate the sodium emission in every frame
; =====================================================================================================================
if part eq 2 or part eq 99 then begin 
;    spectra_cube = MRDFITS(outdir + 'spectra_cube.fits', 0, header, /fscale, /silent )
    img_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header2, /fscale, /silent )
    if part ne 99 then restore, outdir + 'shift_array.sav'           ; contains shift_array, std_devs, aligned_imaging_cube and aligned_spectra_cube
    frame = reform(aligned_spectra_cube[*,*,0])
    s = size(aligned_spectra_cube, /dimensions) 
    si = size(img_cube, /dimensions)
    window, 0, xpos=winpos_x[0],         ypos=winpos_y[0],           xs=s[0], ys=s[1], title='IDL 0 - FRAME'
    window, 1, xpos=winpos_x[0],         ypos=winpos_y[0]+s[1]+40,   xs=s[0], ys=s[1], title='IDL 1 - SCALED REFERENCE'
    window, 2, xpos=winpos_x[0],         ypos=winpos_y[0]+2*s[1]+80, xs=s[0], ys=s[1], title='IDL 2 - EXOSPHERE'
    window, 3, xpos=winpos_x[0]+s[0]+20, ypos=winpos_y[0],           xs=si[0], ys=si[1], title='IDL 3 - FRAME'
;    window, 4, xpos=winpos_x[0]+s[0]+20, ypos=winpos_y[0]+s[1]+40,   xs=s[0], ys=s[1]
    tv, bytscl(Frame)
    
    ;For this data there is no dark or bias calibration frames
    ;Therefore, use the conventional mode date from the following night
    
    Dark = MRDFITS(dark_dir + dark_file, 0, header, /fscale, /silent )
    dark = dark[sps(0):sps(1),sps(2):sps(3)]
    dark = sigma_filter( dark, 5, N_sigma=5)
    acre, dark, dark, thresh, width                 ; try acre, still hot pixels in dark

    if night eq '15' then begin
      spectra_sky_cube = MRDFITS(outdir + 'spectra_sky_cube.fits', 0, header, /fscale, /silent )
      ireference = spectra_sky_cube - dark
      acre, ireference, reference, thresh, width
      tv, bytscl(reference)

      ;get a spectral flat
      iflat = MRDFITS(flat_dir + flat_file, 0, header, /fscale ) 
      iflat = (iflat[sps(0):sps(1),sps(2):sps(3)] - dark)
      iflat = iflat / median(iflat) 
      acre, iflat, flat, thresh, width
;      flat = smart_shift(flat, 0., -.65,  /interp) ;hack
      flat = flat + (1. - median(flat))          ; this normalizes again such that median(flat) = 1
      tv, bytscl(flat)

      reference = reference / flat
      tv, bytscl(reference) 
      
      ; Now find the xoffset in the reference spectrum (likely due to slightly different grating angles)
      y_Mercury = where( total(frame,1) ge 0.5*max(total(frame,1)) )    ; FWHM range of mercury
      frame_spectrum = total(frame(*,y_Mercury),2)
      ref_spectrum = total(reference(*,y_Mercury), 2)
      frame_D1_x = ( where( frame_spectrum(0:s[0]/2) eq min(frame_spectrum(0:s[0]/2))) ) [0]
      ref_D1_x   = ( where( ref_spectrum(0:s[0]/2)   eq min(ref_spectrum(0:s[0]/2))   ) ) [0]   ; hack (we only want the D1 line, so we just check the 1st half of spectrum)

      reference = shift(reference, frame_D1_x - ref_D1_x)
      ref_spectrum = shift(ref_spectrum, frame_D1_x - ref_D1_x)

      exosphere_spectra_cube = fltarr( size(aligned_spectra_cube, /dimensions) )
;stop 
       ;now scale and subtract 
      for i = 0, s[2]-1 do begin                                          ; LOOP over each frame in the series
        frame = reform(aligned_spectra_cube[*,*,i]) - dark                        ;dark subtract the spectrum
        frame = frame / flat                                              ;flatfield the spectrum
      
        ; Now find best scaling for ref_spectrum and subtract
        
        illum = frame / reference                                         ; scale the illumination against Mercury's disk
        for iD = 0, 2, 2 do illum[Na_D_rngs(iD+0):Na_D_rngs(iD+1),*] = !values.F_NaN  ; carefully avoid sodium emissions when scaling to the disk
      
;        illum_along_slit = mean(illum, dimension=1, /nan)
        illum_along_slit = MEDIAN(illum, dimension=1)
        illum_along_slit = MEDSMOOTH( illum_along_slit, 5 )               ;helps with hot pixels from bad flat-fielding where the slit has some dust
        scaled_reference = reference * rebin(transpose(illum_along_slit), s[0], s[1])
        Just_exosphere = Frame - scaled_reference * iref_factor
stop
        ; Check that we've done a good job of minimizing spectrum between Na D lines
        between_D = smooth( just_exosphere(Na_D_rngs(1)+20:Na_D_rngs(2)-20,*), 6, /edge_truncate )
        mean_vals = abs([ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
                          mean(between_D(*,0:min(Y_Mercury))), $              ; mean value below mercury's disk
                          mean(between_D(*,max(Y_Mercury):*)) ])              ; mean value above Mercury's disk
        
        if mean_vals(0) gt 2.*mean(mean_vals(1:2)) then begin             ; find the best "ref_factor" if we haven't
          nfacs = (maxfac - minfac)/dfac + 1                              ; # of factors to consider
          factors = cgscalevector(indgen(nfacs),minfac,maxfac)
          fac_diffs = fltarr(nfacs)                                       ; keep track of the differences in mean_vals for each factor
          for ifac = 0, nfacs-1 do begin
            Just_exosphere = frame - scaled_reference*factors(ifac)
            between_D = smooth( just_exosphere(Na_D_rngs(1)+20:Na_D_rngs(2)-20,*), 6 )
            mean_vals = abs([ mean(between_D(*,min(Y_Mercury):max(Y_Mercury))), $ ; mean value across Mercury's disk
                              mean(between_D(*,0:min(Y_Mercury))), $              ; mean value below mercury's disk
                              mean(between_D(*,max(Y_Mercury):*)) ])              ; mean value above Mercury's disk
            fac_diffs(ifac) = abs( mean_vals(0) - mean(mean_vals(1:2)) )            
          endfor
          best_fac = factors(( where( fac_diffs eq min(fac_diffs) ) )[0])
          print, i, best_fac
          just_exosphere = frame - scaled_reference * best_fac
        endif
        
        wset, 0
        cgimage, bytscl(frame, 0, 200), /axes
        cgoplot, [1,1]*ref_D1_x, !Y.Crange
        cgoplot, [1,1]*frame_D1_x, !Y.Crange, col='red'
        wset, 1
        cgimage, bytscl(scaled_reference, 0, 200), /axes
        cgoplot, [1,1]*ref_D1_x, !Y.Crange
        cgoplot, [1,1]*frame_D1_x, !Y.Crange, col='red'
        wset, 2
        cgimage, bytscl(Just_exosphere, 0, 200), /axes
        cgoplot, [1,1]*ref_D1_x, !Y.Crange
        cgoplot, [1,1]*frame_D1_x, !Y.Crange, col='red'
        wset, 3
        img_frame = shift(reform(img_cube(*,*,i)), [shift_array[i,*]])
        cgimage, bytscl( img_frame, 0, 3500 )
        wait, 0.02
        exosphere_spectra_cube[ *, *, i] = Just_exosphere
        cgtext, 1, 1, 'Frame ' + strfix(1+i)
;        stop
      endfor ;i (loop over frames)
      save, exosphere_spectra_cube, filename = outdir + 'exosphere_spectra_cube.sav'
    endif ;night=='25'

    beep
;   stop
endif

; =====================================================================================================================
; Part 3 : extract and place the 1D spectra.
; =====================================================================================================================
if part eq 3 or part eq 99 then begin 
  if part ne 99 then restore, outdir + 'exosphere_spectra_cube.sav'
  s = size(exosphere_spectra_cube, /dimensions) 
  
  ;Dispersion calculation: at y of 265, lines occur at 390 and 632 in x 
  dispersion = (5895.92424 - 5889.95095) / (632. - 390.) ;Angstroms per pixel
  dispersion_Velocity = 299792.*dispersion / 5889.95095  ;Km/s per pixel at Na D2

   ;-------------------------------------------Extraction----------------------------------------------------------- 

    ;run MPFIT once to find the D2 line center
      img = total(exosphere_spectra_cube, 3)
      search = 16           ;pixels to search over
      expected_pixel = mean(Na_D_rngs(0:1))  ;Rough D2 pixel location
      D2_height = fltarr(n_elements(img[0,*]))
      D2_center = fltarr(n_elements(img[0,*]))
      D2_width  = fltarr(n_elements(img[0,*]))
      for i = 0, n_elements(img[0,*]) - 1 do begin
        result = mpfitpeak(findgen(search*2. + 1), img[expected_pixel-search:expected_pixel+search, i], a, STATUS = STATUS, /positive) 
        if STATUS ge 1 then D2_height[i] = A[0] else D2_height[i] = !values.F_nan
        if STATUS ge 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
        if STATUS ge 1 then D2_width[i] = A[2] else D2_width[i] = !values.F_nan
      endfor
;      stop
      y = findgen(n_elements(img[0,*]))
      height = smooth(D2_height, 5, /edge_mirror) / s[2]
      ;height[0:25] = 1. & height[280:*] = 1.
      height[0:25] = 1. & height[280:*] = 1.
      real = where(finite(D2_center), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
      location = poly(findgen(n_elements(img[0,*])), coeff) + expected_pixel-search 
      real = where(finite(D2_width), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D2_width[real], 2)
      width = poly(findgen(n_elements(img[0,*])), coeff) 
      dummy = img
      dummy[location, y] = max(dummy) 
      tv, bytscl(dummy, -100, 200)

    ;setup MPFIT Gaussian parameters
      guess = [80.,location[0],1.5,0.]
      A = guess
      parinfo = replicate( {value:0.D, fixed: 0b, limited: [0b,0b], limits: dblarr(2) }, 3)
      parinfo[0].limited = 1b                                 ;limit amplitude 
      parinfo[0].limits = [0.1, 1000. ]                       ;positive only
      parinfo[1].limited = 1b                                 ;limit sigma   
      parinfo[1].limits = [14., 18.]                          ;limit sigma width in pixels      
      parinfo[2].limited = 1b                                 ;limit sigma   
      parinfo[2].limits = [1., 5.]                            ;limit sigma width in pixels
  
  brightness = fltarr(s[2], s[1]) & err_brightness = fltarr(s[2], s[1])
  linewidth = fltarr(s[2], s[1]) & err_linewidth = fltarr(s[2], s[1])
  for n = 0, s[2]-1 do begin
;      img = exosphere_spectra_cube[*,*,n]
      img = smooth( exosphere_spectra_cube[*,*,n], [2,6] )    ; hack (not in AEOS version, we want smoother spatial structure in spectrum)
      err_img = sqrt( abs(img) )     
      for i = 0, s[1]-1 do begin
               a[0] = height[i]
               a[1] = location[i] - expected_pixel + search
               a[2] = Width[i]
               
;               cgplot, findgen(search*2. + 1), img[location[i]-search:location[i]+search, i]
               result = mpfitpeak(findgen(search*2. + 1), img[location[i]-search:location[i]+search, i], A, PERROR = Err_A, $
                                  error = err_img[location[i]-search:location[i]+search, i], /positive, /NAN, STATUS = STATUS, PARINFO = PARINFO, NTERMS =3)
;               cgoplot, findgen(search*2. + 1), result, color = 'red'
;               stop
               if STATUS ge 1 then begin
                  brightness[n,i] = A[0]*A[2]*SQRT(2*!DPI) 
                  err_brightness[n,i] = brightness[n,i] * sqrt( (err_A[0]/A[0])^2. + (err_A[2]/A[2])^2. )
                  linewidth[n,i] = dispersion_Velocity*sqrt((2.355*A[2])^2 - (2.45)^2) ;FWHM KM/S instrumental measured width is 2.45 pixels 
                  err_linewidth[n,i] = dispersion_Velocity*2.355*Err_A[2]
               endif else begin
                  brightness[n,i] = !values.F_nan
                  linewidth[n,i] = !values.F_nan
               endelse 
      endfor
  endfor  
  tv, bytscl(brightness, 0, 1000)
  save, brightness, linewidth, filename = outdir + 'brightness.sav'
  beep
;  stop
endif  

; =====================================================================================================================
; Part 4 : put everything together and build images and movies of the exosphere.
; =====================================================================================================================
if part eq 4 or part eq 99 then begin  
   if part ne 99 then restore, outdir + 'brightness.sav'            ; contains brightness, linewidth
   if part ne 99 then restore, outdir + 'shift_array.sav'           ; contains shift_array, std_devs, and aligned_imaging_cube
   if part ne 99 then imaging_cube = MRDFITS(outdir + 'imaging_cube.fits', 0, header, /fscale, /silent )
   s = size(imaging_cube, /dimensions)
   
   brightness = congrid(brightness, s[2], s[1]) ;HACK NEED PROPER SCALING OF THE SPECTRUM TO IMAGING PLATESCALES
   home = [341,238] ;Mercury centroid of the last frame number 500, to which all the other frames are aligned
   home = [236,165]
   home = [231,165]
;   home = [131,102]  ; trying home values for a smaller cutout of the imaging and spectral frames to speed things up
   
;   big_cube = fltarr(s[0], 2*s[1], s[2])
   big_cube_spec = fltarr(s[0], s[1], s[2])
;   big_cube(*,*,*) = !Values.F_NAN
   big_cube_spec(*,*,*) = !Values.F_NAN
   img_cube = fltarr(s[0], s[1])
   img_count = fltarr(s[0], s[1])                ; array to track number of "real" pixels (ie not obscured by slit)
   spec_cube = fltarr(s[0], s[1])
   spec_dwell = fltarr(s[0], s[1])                                        ; keep track of the # of frames in each pixel
   spec_temp  = fltarr(s[0], s[1]) 
   big_frame_total = fltarr(s[0], 2*s[1])
   window, 0, xs = s[0], ys = 2*s[1], xpos=winpos_x[0], ypos=winpos_y[0]
;   cgdisplay, s[0], s[1]*2, xpos=winpos_x[0], ypos=winpos_y[0], wid=0
;   window, 1, xs = s[0], ys = s[1], xpos=winpos_x[2], ypos=winpos_y[2]
   xs_spec = spec_extraction(1) - spec_extraction(0) + 1
   ys_spec = spec_extraction(3) - spec_extraction(2) + 1
   xs_img= img_extraction(1) - img_extraction(0) + 1
   ys_img = img_extraction(3) - img_extraction(2) + 1
   window, 1, xs = xs_img*movie_scale, ys = ys_img*movie_scale, xpos=winpos_x[0]+s[0]+20, ypos=winpos_y[0]   
;   cgdisplay, xs_spec*movie_scale, ys_spec*movie_scale, xpos=winpos_x[0]+40, ypos=winpos_y[0], wid=1
;   window, 2, xs = s[0], ys = s[1], xpos=winpos_x[2], ypos=winpos_y[2]+s[1]+35
   window, 2, xs = xs_img*movie_scale, ys = ys_img*movie_scale, xpos=winpos_x[0]+s[0]+20, ypos=winpos_y[0]+ys_img*movie_scale+40
   window, 3, xs = s[0], ys = 2*s[1], xpos=winpos_x[0], ypos=winpos_y[0]+2*s[1]+40
;stop
   ; Some video related variables
   mpgFilename = outdir + 'Mercury Na image.mp4'
   mpgFilename2 = outdir + 'Mercury disk.mp4'
   mpgFilename3 = outdir + 'Mercury scan.mp4'
   video = IDLffVideoWrite(mpgFilename, format='mp4')
   video2 = IDLffVideoWrite(mpgFilename2, format='mp4')
   video3 = IDLffVideoWrite(mpgFilename3, format='mp4')
   framerate = 30
   framedims = [xs_img*movie_scale, ys_img*movie_scale]
   framedims2 = [xs_img*movie_scale, ys_img*movie_scale]
   framedims3 = [s[0], 2*s[1]]
   stream  = video.AddVideoStream(framedims[0], framedims[1], framerate)
   stream2 = video2.AddVideoStream(framedims2[0], framedims2[1], framerate)
   stream3 = video3.AddVideoStream(framedims3[0], framedims3[1], framerate)

   good_frames = where(std_devs gt stddev_cutoff, complement=bad_frames)  ; these are the good frames we'll use below (in terms of good seeing)
   ngood_shifts = n_elements(where(shift_array(*,0) ne -666))             ; these are the good frames (in terms of a reasonable shift_array value)

;   for ig = 253, n_elements(good_frames)-1 do begin                         ; loop over good frames
   for ig = 0, n_elements(good_frames)-1 do begin                         ; loop over good frames
    
        i = good_frames(ig)
        if shift_array[i,0] eq -666 then goto, skipbad
        
;   for i = 0, 20 do begin;s[2]-1 do begin
        big_frame = fltarr(s[0], 2*s[1])
        empty_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        empty_frame[home[0]+shift_array[i,0]-1:home[0]+shift_array[i,0]+1,*] = [brightness[i, *], brightness[i, *], brightness[i, *]] ;Make it bigger?
        big_frame[*, s[1]:*] = rotate(empty_frame, 7) ;* 5. ;the factor of 5 is a HACK FOR NOW
        big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i])
        big_frame_total = big_frame_total + big_frame
;        stop
        wset, 0
        cgloadct, ct, /silent
        big_frame_scaled = big_frame / (total(big_frame,/nan)/7.427e7)
        cgimage, bytscl(big_Frame_scaled, 0, 10000)   ; 7e7 is hack for total of all frames
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_scale
        cgtext, 0.02, 0.96, 'Extracted Na emission', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='black', /normal, charthick=1.25*movie_Scale
        cgtext, 0.02, 0.45, 'IMAGING CHANNEL', charsize=0.5*movie_Scale, col='white', /normal, charthick=0.5*movie_Scale
        cgtext, .75, .97, 'Frame # ' + strfix(i+1), /data, col='black', /normal, charthick=1.5
        movie_png = cgsnapshot(filename = outdir + 'Mercury scan.png', /nodialog)
        image = read_png(outdir + 'Mercury scan.png')
        void = video3 -> Put(stream3, image)
               
;        big_cube[*,*,i] = big_frame
        spec_cut_for_big_cube_spec = reform(big_frame(*,s[1]:*))
        big_cube_spec[*,*,i] = smooth( spec_cut_for_big_cube_spec, [0,8], /nan )
        img_cube = img_cube + reform(big_frame(*,0:s[1]))
        spec_frame = reform(big_frame(*,s[1]:*))
        spec_real = where(finite(spec_frame),complement=spec_bkg)         ; find elements with real values, otherwise assign to background
        spec_summed = total(finite(spec_frame),2)                         ; by summing rows we can find the columns corresponding to the slit location
        slit_columns = where(spec_summed eq max(spec_summed))             ; columns corresponding to slit location
        if n_elements(slit_columns) le 3 then begin                       ; sometimes we have a bad frame, in which case slit_columns thinks everything is beneath the slit
          img_count = img_count + 1.                                        ; increment "real" pixels
          img_count(slit_columns,*) = img_count(slit_columns,*) - 1.        ; but also account for the pixels covered by the slit
;        spec_dwell(spec_real) = spec_dwell(spec_real) + integration       ; dwell time
          spec_dwell(min(slit_columns):max(slit_columns),*) = spec_dwell(min(slit_columns):max(slit_columns),*) + integration
        endif
        spec_frame(spec_bkg) = 0d
        spec_frame = smooth(spec_frame, [0,8])                            ; HACK (just in the sense that maybe we don't want to smooth here)
        spec_cube(spec_real) = spec_cube(spec_real) + spec_frame(spec_real)

;wset, 1
;cgimage, spec_frame
;wset, 2
;cgimage, spec_dwell
;stop

        wset, 1
        spec_dwell_num = where(spec_dwell gt 0.)                          ; find pixels with spectral information included
        ;if i lt s[2]-1 then $
        if ig lt n_elements(good_frames)-1 then $
          spec_temp(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num) + spec_frame(spec_dwell_num) else $   
          spec_temp(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num)  ; ie just don't plot slit contribution for last frame
        spec_frame_slit = spec_frame(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
        spec_frame_slit_value = where(spec_frame_slit gt 0., complement=spec_frame_slit_bkg)
        spec_frame_slit_2 = spec_frame_slit
        spec_frame_slit_2(spec_frame_slit_bkg) = 0.;!Values.F_NAN
        blah = big_cube_spec(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3),*)
;        blah = big_cube_spec(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3),*)
;        wslit = where( total(spec_frame_slit,2) eq max(total(spec_frame_slit,2)) )
        if do_median then spec_movie_temp = median(blah, dimension=3, /even) $
                     else spec_movie_temp = mean(blah, dimension=3, /nan)
        spec_movie = bytscl(spec_movie_temp, max_spec*0.05, max_spec*.8)
;        spec_movie = bytscl(spec_temp(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3)), max_spec*0.05, max_spec*.8)
;        spec_movie = congrid(spec_movie, xs_spec*movie_scale, ys_spec*movie_scale)
;        if n_elements(wslit) eq 3 then spec_movie(min(wslit)>1:max(wslit)<(n_elements(spec_movie(*,0))-1),*) = 1e3 ; hack to make slit show up
        axis_format = {XTicklen:0.0001, YTickLen:0.0001, Xthick:4, Ythick:4}
        if do_rotate then cgimage, rotate(spec_movie,5), /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format $
                     else cgimage, spec_movie, /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format
        if do_rotate then cgimage, rotate(spec_frame_slit_2,5), transparent=20, missing_Value=0.0, ctindex=0 $
                     else cgimage, spec_frame_slit_2,  transparent=20, missing_Value=0.0, ctindex=0
        img_cube_cut = img_cube(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
        cgcontour, normalize(img_cube_cut), /onimage, levels=[0.2], label=0

;if ig gt 31 then stop        
        ; Next 2 lines draw slit
;        cgoplot, [1,1]*min(wslit), !Y.Crange, line=0, thick=2
;        cgoplot, [1,1]*max(wslit), !Y.Crange, line=0, thick=2

;        cgcolorbar, /vertical, charsize=1.e-3, position=[0.09,0.1,0.13,0.8];, oob_low='white'
        cgcolorbar, /vertical, charsize=1.e-3, position=[0.9,0.15,0.95,0.85]
        cgtext, .65, .05, '# frames = ' + strfix(ig+1), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
        cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
        movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
        image = read_png(outdir + 'Mercury Na.png')
        void = video -> Put(stream, image)
;stop
        wset, 2
;        img_movie = bytscl(img_cube(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3)), max_img*0.2, max_img)
        img_count_cut = smooth(img_count(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3)),6,/edge_truncate)
;        img_movie = bytscl(img_cube_cut, 0.2*max(img_cube_cut), max(img_cube_cut))
;        img_movie = congrid(img_movie, xs_img*movie_scale, ys_img*movie_scale)
        img_movie = congrid(img_cube_cut * (count_factor*(max(img_count_cut)/img_count_cut)^2.), xs_img*movie_scale, ys_img*movie_scale)
        cgloadct, 0, /silent
        cgimage, img_movie
        cgcolorbar, /vertical, charsize=1.e-3
        cgtext, .65, .05, '# frames = ' + strfix(ig+1), /data, /normal, charthick=0.375*movie_Scale, charsize=0.5*movie_scale
        cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale
        movie_png = cgsnapshot(filename = outdir + 'Mercury disk.png', /nodialog)
        image = read_png(outdir + 'Mercury disk.png')
        void = video2 -> Put(stream2, image)

        wset, 3
        cgimage, bytscl(big_frame_total/double(i),0,5000)
        skipbad: wait, 0.02   

;loadct, 20
;cgimage, reform(aligned_imaging_cube(*,*,i)), /axes, /keep_aspect_ratio
;x = spec_dwell*0.
;x(slit_columns,*) = 1.0
;cgcontour, x, /onimage, label=0
;print, slit_columns
;stop
   endfor ;ig
   
   ; LOOP over bad frames too for comparison figure ???
   dobad = 1
   if dobad then begin
     img_cube_bad = fltarr(s[0], s[1])
     for ib = 0, n_elements(bad_frames)-1 do begin                         ; loop over bad frames
        i = bad_frames(ib)
        if shift_array[i,0] eq -666 then goto, skipbad2

        big_frame = fltarr(s[0], 2*s[1])
        empty_frame = make_array(s[0], s[1], /float, value = !values.F_NaN)
        empty_frame[home[0]+shift_array[i,0]-1:home[0]+shift_array[i,0]+1,*] = [brightness[i, *], brightness[i, *], brightness[i, *]] ;Make it bigger?
        big_frame[*, s[1]:*] = rotate(empty_frame, 7) ;* 5. ;the factor of 5 is a HACK FOR NOW
        big_frame[*, 0:s[1]-1] = reform(aligned_imaging_cube[*,*,i < ngood_shifts-1])
        big_frame_total = big_frame_total + big_frame
;        big_cube[*,*,i] = big_frame
        img_cube_bad = img_cube_bad + reform(big_frame(*,0:s[1]))

        skipbad2: wait, 0.02
      endfor ;ib     
    img_cube_cut_bad = img_cube_bad(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
   endif ;dobad==1
     
   
   wset, 3
;   test = median(big_cube, dimension = 3, /even)
   ;tv, bytscl(smooth(test,[3,3], /NaN), 0, 10000)
;   cgimage, bytscl(test, 0, 4500)
   spec_dwell_num = where(spec_dwell gt 0.)                               ; find pixels with spectral information included
   spec_cube(spec_dwell_num) = spec_cube(spec_dwell_num) / spec_dwell(spec_dwell_num)  ; normalize by # of frames summed
   
   ; Now add some blank frames to the Na image movie before plotting a smoothed frame
   ;spec_raw = spec_cube(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
   wset, 1
   spec_raw = spec_cube(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
   spec_raw = congrid(spec_raw, xs_spec*movie_scale, ys_spec*movie_scale) ; embiggen
;   spec_movie = bytscl(spec_raw, max_spec*0.05, max_spec*.8)
   loadct, ct, /silent
   if do_rotate then cgimage, rotate(spec_movie,5), /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format $
                else cgimage, spec_movie, /axes, charsize=1e-3, position=[0.0,0.0,1.,1.], axkeywords=axis_format
   cgcontour, normalize(img_cube_cut), /onimage, levels=[0.2], label=0
   cgcolorbar, /vertical, charsize=1.e-3, position=[0.9,0.15,0.95,0.85]
   cgtext, .65, .05, '# frames = ' + strfix(ig), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
   cgtext, 0.02, 0.96, 'Mercury Na', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
   movie_png = cgsnapshot(filename = outdir + 'Mercury Na.png', /nodialog)
   image = read_png(outdir + 'Mercury Na.png')

   for i = 0, 30 do void = video -> Put(stream, image)
   smooth_width = 4
   if do_rotate then cgimage, rotate(smooth(spec_raw,smooth_width*movie_scale,/edge_truncate), 5) $              ; HACK
                else cgimage, smooth(spec_raw,smooth_width*movie_scale,/edge_truncate)
   cgcontour, normalize(img_cube_cut), /onimage, levels=[0.2], label=0
   cgtext, .65, .05, '# frames = ' + strfix(ig), /normal, charthick=0.5*movie_Scale, font=-1, charsize=0.5*movie_Scale
   cgtext, 0.02, 0.96, 'Mercury Na (smoothed=' + strfix(smooth_width) + ')', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale, font=-1
   movie_png = cgsnapshot(filename = outdir + 'Mercury Na - smoothed.png', /nodialog)
   image = read_png(outdir + 'Mercury Na - smoothed.png')
   for i = 0, 30 do void = video -> Put(stream, image)
      
   video -> Cleanup
   video2 -> Cleanup
   video3 -> Cleanup

   ; Overplot Na contours on top of Mercury disk
   loadct, 0, /silent
   cgimage, bytscl(img_movie), /axes
   spec_smoothed = bytscl(smooth(spec_raw,smooth_width*movie_scale,/edge_truncate))
;   cgcontour, spec_smoothed, /onimage, color='black', levels=[50,100,150,175,200,225,250], xthick=3 
;   cgcontour, spec_smoothed, /onimage, color='gold', levels=[50,100,150,175,200,225,250]
   cgcontour, normalize(double(spec_smoothed)), /onimage, color='gold', levels=[0.2,0.4,0.6,0.8]
   al_legend, ['Disk','Na'], box=0, /bottom, /right, charsize=0.5*movie_scale, line=[0,0], linsize=0.2, thick=[6,2], color=['gray','gold'], textcolor='white'
   combined = cgsnapshot(filename = outdir + 'Mercury Na + disk.png', /nodialog)

   ; Generate the "bad images" Mercury disk for comparison
   if dobad then begin
     bad_img = bytscl(img_cube_cut_bad, 0.2*max(img_cube_cut_bad), max(img_cube_cut_bad))
     cgimage, bad_img
;    cgimage, img_cube_cut_bad
     cgcolorbar, /vertical, charsize=1.e-3
     cgtext, .65, .05, '# frames = ' + strfix(ib), /data, /normal, charthick=0.375*movie_Scale, charsize=0.5*movie_Scale
     cgtext, 0.02, 0.96, 'Solar continuum', /normal, charthick=0.75*movie_Scale, charsize=0.625*movie_Scale
     bad_png = cgsnapshot(filename = outdir + 'Mercury disk - bad frames.png', /nodialog)

     ; And might as well compare the good vs bad contours while we're here...
     cgimage, img_movie, /axes
;     cgcontour, img_movie, /onimage, color='gold'
     cgcontour, normalize(img_movie), /onimage, levels=[0.2,0.4,0.6,0.8], color='gold'
;     cgcontour, bad_img, /onimage, color='red'
     cgcontour, normalize(img_cube_cut_bad), /onimage, levels=[0.2,0.4,0.6,0.8], color='red'
     al_legend, ['good','bad'], box=0, /bottom, /right, charsize=2, line=[0,0], linsize=0.2, thick=[4,4], color=['orange','red'], textcolor=['orange','red']
     compare = cgsnapshot(filename = outdir + 'Mercury disk - good vs bad contours.png', /nodialog)
   endif ;dobad==1
   

   
   ; Create image showing the spatial coverage of the slit position
   wset, 2
;   dwell_img = spec_dwell(spec_extraction(0):spec_extraction(1),spec_extraction(2):spec_extraction(3))
   dwell_img = spec_dwell(img_extraction(0):img_extraction(1),img_extraction(2):img_extraction(3))
   loadct, 22, /silent
;   cgloadct, 22, ncolors=max(spec_dwell)/integration, bottom=0, /silent
   cgimage, bytscl(dwell_img, 0., max(dwell_img)), /axes, title='Slit Coverage'
   cgcontour, normalize(img_cube_cut), /onimage, levels=[0.2]
   cgloadct, 22, ncolors=max(spec_dwell)/integration, bottom=0, /silent  
   cgcolorbar, /vertical, maxrange=max(spec_dwell), title='Total Na data coverage (sec)', /discrete, ncolors=max(spec_dwell)/integration, position=[0.9,0.2,0.95,0.8], tickinterval=integration*4.
   combined = cgsnapshot(filename = outdir + 'Mercury slit coverage.png', /nodialog)
   
   
   
   beep
   stop
endif

; cgPS2Raster, 'test_2.ps', /JPEG
;   IDL> cgPS_Open, 'test_2.ps'
;   IDL> Graphic_Display
;   IDL> cgPS_Close, /PNG
  
end



