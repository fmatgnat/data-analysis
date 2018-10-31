pro RIPS_Mercury_Perkins_reference

; *********************************************************************************************************************
; *********************************************************************************************************************
; Routine to create a reference Mercury image that RIPS_Mercury_Perkins_v03.pro can use for image alignment.
;
; Input variables
; 
; Created:  29 Oct 2018 - 
; Edits:        
; *********************************************************************************************************************
; =====================================================================================================================
; Define variables
; =====================================================================================================================
input_dir     = 'Y:\obs_18\Perkins_RIPS_March\15\renamed\'
input_series  = 'Mercury_Na_804.fits'
dark_dir      = 'Y:\obs_18\Perkins_RIPS_March\14\Carl_Keep\'              ; directory with dark file
flat_dir      = 'Y:\obs_18\Perkins_RIPS_March\15\'                        ; directory with flat file
dark_file     = 'Dark.fits'                                               ; dark frame
flat_file     = 'Flat_NaSpectra_624slitwidth_NaND1imaging1.fits'          ; flat frame
width         = 12
thresh        = 5
smooth_width  = 6
scale         = 1
frame_max     = 1e6
ims           = [283,759,39,374]                                          ; ranges for image extraction [x1,x2,y1,y2]

DEVICE, GET_SCREEN_SIZE = screen_size                                     ; find screen size
winpos_x      = screen_size(0)/2                                          ; default initial x window location
winpos_y      = 0                                                         ; default initial y window location
!p.multi      = 0.
!p.position   = 0.
cgloadct, 20, /silent

; =====================================================================================================================
; Read in, co-align and save master reference
; =====================================================================================================================
icube         = readfits( input_dir + input_series, header, /silent)
iflat         = readfits( flat_dir + flat_file, flat_header, /silent)
idark         = readfits( dark_dir + dark_file, dark_header, /silent)
s             = size(icube,/dimension)
iint          = sxpar(header, 'EXPOSURE')                                 ; integration time for individual frames
dint          = sxpar(dark_header, 'EXPOSURE')                            ; dark integration time
fint          = sxpar(flat_header, 'EXPOSURE')                            ; flat integration time

for i = 0, s[2]-1 do begin

  if i eq 0 then begin
    image     = reform(icube(*,0:s[1]/2-1,i)) / iint
    x_center  = ( where( total(image,2) eq max(total(image,2)) ) )[0]
    y_center  = ( where( total(image,1) eq max(total(image,1)) ) )[0]
    x_range   = where( normalize(total(image,2) - mean(total(image,2))) gt 0.03 )
    y_range   = where( normalize(total(image,1) - mean(total(image,1))) gt 0.03 )
;    ims       = [min(x_range)-20,max(x_range)+20,min(y_range)-20,max(y_range)+20]
    dark      = reform(idark[ims(0):ims(1),ims(2):ims(3),0]) / dint       ; divide by dint to normalize to dark counts/sec
    iflat2    = iflat[ims(0):ims(1),ims(2):ims(3)]                        ; extract flat 
;    flat      = ( iflat2 - dark ) / max( iflat2 - dark )
    flat      = iflat2 / max(iflat2)
    flat      = flat + (1. - median(flat))                                ; this normalizes such that the median of the flat is 1
    iref      = image[ims(0):ims(1),ims(2):ims(3)] / iint                 ; extract frame and normalize to counts/sec
    reference = (iref - dark) * flat                                      ; correct for flat field structure
    timg      = reference * 0.
    window, 0, xs = (ims[1]-ims[0])*scale, ys = (ims[3]-ims[2])*scale, xpos=winpos_x, ypos=winpos_y
    window, 1, xs = (ims[1]-ims[0])*scale, ys = (ims[3]-ims[2])*scale, xpos=winpos_x, ypos=winpos_y+(ims[3]-ims[2])*scale+40
    window, 2, xs = (ims[1]-ims[0])*scale, ys = (ims[3]-ims[2])*scale, xpos=winpos_x+(ims[1]-ims[0])*scale+20, ypos=winpos_y
    window, 3, xs = (ims[1]-ims[0])*scale, ys = (ims[3]-ims[2])*scale, xpos=winpos_x+(ims[1]-ims[0])*scale+20, ypos=winpos_y+(ims[3]-ims[2])*scale+40
  endif else begin
    image     = reform(icube(*,0:s[1]/2-1,i)) / iint
    iframe    = image[ims(0):ims(1),ims(2):ims(3)] / iint                 ; extract frame and normalize to counts/sec
    frame     = (iframe - dark) * flat                                    ; correct for flat field structure
    CORREL_OPTIMIZE, reference, frame, xoff, yoff, NUMPIX=150
    aframe    = shift(frame, [xoff,yoff])
    x_center2 = ( where( total(reference,2) eq max(total(reference,2)) ) )[0]
    y_center2 = ( where( total(reference,1) eq max(total(reference,1)) ) )[0]
    if abs(xoff) lt 50 and abs(yoff) lt 50 then begin                     ; ie don't add the frame if the offset failed
      print, 'Frame ' + strfix(i) + ' offsets (x,y):', xoff, yoff
      timg = timg + aframe
    endif else print, 'CORREL failed for frame ' + strfix(i)
    
    ; And plot results
    wset, 0
    cgimage, bytscl(iframe,0,frame_max), /axes, /keep_aspect, title='raw frame'
    cgcontour, normalize(smooth(iframe,smooth_width)), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='black'
    cgoplot, !X.Crange, [1,1]*y_center2, col='blue'
    cgoplot, [1,1]*x_center2, !Y.Crange, col='blue'
    wset, 1
    cgimage, bytscl(frame,0,frame_max), /axes, /keep_aspect, title='dark+flat corrected frame'
    cgcontour, normalize(smooth(frame,smooth_width)), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='black'
    cgoplot, !X.Crange, [1,1]*y_center2, col='blue'
    cgoplot, [1,1]*x_center2, !Y.Crange, col='blue'
    wset, 2
    cgimage, bytscl(reference,0,frame_max), /axes, /keep_aspect, title='reference'
    cgcontour, normalize(smooth(reference,smooth_width)), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='black'
    cgoplot, !X.Crange, [1,1]*y_center2, col='blue'
    cgoplot, [1,1]*x_center2, !Y.Crange, col='blue'
    wset, 3
    cgimage, bytscl(aframe,0,frame_max), /axes, /keep_aspect, title='aligned frame: ' + strfix(i)
    cgcontour, normalize(smooth(reference,smooth_width)), /onimage, levels=[0.2,0.4,0.6,0.8], c_color='black'
    cgoplot, !X.Crange, [1,1]*y_center2, col='blue'
    cgoplot, [1,1]*x_center2, !Y.Crange, col='blue'
;    stop
  endelse
  
  wait, 0.05
endfor ;i

output_file   = ( strsplit(input_series,'.',/extract) )[0] + '_reference_image.fits'
MWRFITS, timg, input_dir + output_file, header, /create, /silent 

stop

end