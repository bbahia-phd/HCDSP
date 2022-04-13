# Function RotSurvey
#
# Purpose: Rotate the given coordinates to horizontal
#
#
#  Inputs:
#        vx and vy - vectors with x and y coordinates
#
#  Outputs:
#        vxr and vyr - vectors with rotated x and y coordinates
#
# Breno Bahia hopes this is correct.
##

function rotate_survey(vx,vy)

   # Get centre-of-mass
   vxc = mean(vx);
   vyc = mean(vy);
   
   # Translate survey so COM is new origin
   vxx = vx .- vxc;
   vyy = vy .- vyc;
   
   # Find angle to rotate
   ang = (vy[end]-vy[1])/(vx[end]-vx[1]);
   ang = -atan(ang)

   # Rotation matrix brute forte
   c = cos(ang);  
   s = sin(ang);

   # Actual rotation
   vxr = vxx.*c .- vyy.*s;
   vyr = vxx.*s .+ vyy.*c;
   
   # Translate back
   vxr = vxr .+ vxc;
   vyr = vyr .+ vyc;   

   # Return rotated vectors
   return vxr,vyr,ang

end
