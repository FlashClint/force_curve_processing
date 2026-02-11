# force_curve_processing
a MATLAB script to batch-process AFM approach force curves

Force curve data format (.xlsx):
- Approach portion only, baseline flattened, calibrated
- Excel layout: columns alternate: col1 = distance (nm), col2 = force (nN), col3 = distance (nm), col4 = force (nN), etc.

Model:
- Hertz spherical indenter: F = (4/3) * E/(1-nu^2) * sqrt(R) * delta^(3/2)
- Bottom-effect (finite-thickness) correction applied using Garcia/Chiodini et al.
   psi(chi) = 1 + 1.133*chi + 1.283*chi^2 + 0.769*chi^3 + 0.0975*chi^4
   where chi = a/h and a = sqrt(R*delta) is the contact radius.

Interactive workflow (per curve):
 1) Show one force-vs-distance curve (force in nN on y, distance in nm on
    x) at a time
 2) User drags two vertical lines to define the fitting region.
 3) Script fits E (Young's modulus) using corrected Hertz model on selected region
    and overlays fitted curve and reports E on the plot.
 4) User can choose to Save this curve and move to next curve or reselect/exit.
