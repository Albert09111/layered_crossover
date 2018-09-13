addpath('lib');

base_info = loadCSVabsFileRawtoStructure('design_file/base_strand_info.csv','');

crossover_info = loadCSVabsFileRawtoStructure('design_file/crossover_info.csv','');
crossover_info(1,:) = [];

scaffold_helix1 = {};
scaffold_helix2 = {};

helper_helix1 = {};
helper_helix2 = {};

for c1 = 1:6
    scaffold_1 = crossover_info{c1,2};
    scaffold_2 = crossover_info{c1,3};
    
    helper_1 = crossover_info{c1,4};
    helper_2 = crossover_info{c1,5};
    
    scaffold_1_pos = get_position(scaffold_1);
    scaffold_2_pos = get_position(scaffold_2);
    
    scaffold_lx = (scaffold_1_pos+scaffold_2_pos)/2;
    
    helper_1_pos = get_position(helper_1);
    helper_2_pos = get_position(helper_2);
    
    helper_lx = (helper_1_pos+helper_2_pos)/2;
    
    scaffold_helix1{end+1} = scaffold_lx;
    helper_helix1{end+1} = helper_lx;
end

for c1 = 7:12
    scaffold_1 = crossover_info{c1,2};
    scaffold_2 = crossover_info{c1,3};
    
    helper_1 = crossover_info{c1,4};
    helper_2 = crossover_info{c1,5};
    
    scaffold_1_pos = get_position(scaffold_1);
    scaffold_2_pos = get_position(scaffold_2);
    
    scaffold_lx = (scaffold_1_pos+scaffold_2_pos)/2;
    
    helper_1_pos = get_position(helper_1);
    helper_2_pos = get_position(helper_2);
    
    helper_lx = (helper_1_pos+helper_2_pos)/2;
    
    scaffold_helix2{end+1} = scaffold_lx;
    helper_helix2{end+1} = helper_lx;
end

%%%% scaffold calculation
scaffold_angle = zeros(6,6);

for c2 = 1:6
    
    horizental_vect = scaffold_helix1{1}-scaffold_helix1{3};
    scaffold_lx1_pos = scaffold_helix1{c2};
    
    for c3 = 1:6
        scaffold_lx2_pos = scaffold_helix2{c3};

        scaffold_vect = scaffold_lx1_pos-scaffold_lx2_pos;
        angle = 2*angle_cal(horizental_vect,scaffold_vect);
        
        if angle <90
            angle = angle;
        elseif angle >90
            angle = 180-90;
        end
            
        
        scaffold_angle(c2,c3) = angle;
        %disp(angle);
    end
    
end


%%% helper calculation
helper_angle = zeros(6,6);
for c2 = 1:6
    horizental_vect = helper_helix1{1}-helper_helix1{3};
    scaffold_help1_pos = helper_helix1{c2};
    for c3 = 1:6
    scaffold_help2_pos = helper_helix2{c3};
    
    scaffold_vect = scaffold_help1_pos-scaffold_help2_pos;
    angle = 2*angle_cal(horizental_vect,scaffold_vect);
    
    if angle <90
        angle = angle;
    elseif angle >90
        angle = 180-90;
    end
    
    helper_angle(c2,c3) = angle;
    %disp(angle);
    
    end
end


%%% hybrid calculation
hybrid_angle = zeros(6,6);

for c2 = 1:6
    
    horizental_vect = scaffold_helix1{1}-scaffold_helix1{3};
    scaffold_lx1_pos = scaffold_helix1{c2};
    
    for c3 = 1:6
        scaffold_lx2_pos = helper_helix2{c3};

        scaffold_vect = scaffold_lx1_pos-scaffold_lx2_pos;
        angle = 2*angle_cal(horizental_vect,scaffold_vect);
        
        if angle <90
            angle = angle;
        elseif angle >90
            angle = 180-90;
        end
            
        
        hybrid_angle(c2,c3) = angle;
        %disp(angle);
    end
    
end



Xtile_lable = {'1','2','3','4','5','6'};
Ytile_lable = {'1','2','3','4','5','6'};

mkdir('figure_output');
close all;
figure(1);
imagesc(scaffold_angle);
colorbar;
caxis([5 95])
axis square;
set(gca,'fontsize',24);
title('strainght strand','fontsize',28);
print(1,'-djpeg','-r600','figure_output/straing_strand_design_space.jpeg')

figure(2);
imagesc(helper_angle);
colorbar;
caxis([5 95]);
axis square;
set(gca,'fontsize',24);
title('crossover strand','fontsize',28);
print(2,'-djpeg','-r600','figure_output/crossover_strand_design_space.jpeg')

figure(3);
imagesc(hybrid_angle);
colorbar;
caxis([5 95]);
axis square;
set(gca,'fontsize',24);
title('crossovers on two types of strand','fontsize',28);
print(3,'-djpeg','-r600','figure_output/hybrid_strand_design_space.jpeg')


















