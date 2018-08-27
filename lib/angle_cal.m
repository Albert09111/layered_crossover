function angle = angle_cal(v1,v2)
        v1_norm = v1/sqrt(v1*v1');
        v2_norm = v2/sqrt(v2*v2');
        
        multiplicatoin = sum(v1_norm.*v2_norm);
        
        angle = acos(multiplicatoin)/pi*180;
        
        if angle < 90
            angle = angle;
        else
            angle = 180-angle;
        end
end