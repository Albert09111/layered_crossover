FileName = 'design_file/double_crossover_tile.dnajson';
fid = fopen(FileName,'r');
data_set_counter = 0;
temp_string = fgetl(fid);
base_info = {'id','postionx','positiony','positionz','basetype','across','up','down'};
base_set = {'Cytosine','Guanine','Thymine','Adenine'};

n = 2;

while temp_string ~= -1;
    
    temp_string = fgetl(fid);

    if ~isempty(strfind(temp_string,'id'))
        id_index = regexp(temp_string,'\d');
        id = temp_string(id_index(1):id_index(end));
        base_info{n,1} = id;
        
    elseif ~isempty(strfind(temp_string,'position'))
        position_index1 = regexp(temp_string,'[');
        position_index2 = regexp(temp_string,',');
        position_index3 = regexp(temp_string,']');
        
        position_x = temp_string(position_index1+1:position_index2(1)-1);
        position_y = temp_string(position_index2(1)+1:position_index2(2)-1);
        position_z = temp_string(position_index2(2)+1:position_index3-1);
        base_info{n,2} = position_x;
        base_info{n,3} = position_y;
        base_info{n,4} = position_z;
        
    elseif ~isempty(strfind(temp_string,'type'))
        if strfind(temp_string,base_set{1});
            type = 'Cytosine';
        elseif strfind(temp_string,base_set{2})
            type = 'Guanine';
        elseif strfind(temp_string,base_set{3})
            type = 'Thymine';
        elseif strfind(temp_string,base_set{4})
            type = 'Adenine';
        end
        base_info{n,5} = type;
        
    elseif ~isempty(strfind(temp_string,'across'))
        across_index = regexp(temp_string,'\d');
        is_null = regexp(temp_string,'null');
        if ~isempty(is_null)
            across_id = 'null';
        else
            across_id = temp_string(across_index(1):across_index(end));
        end
        base_info{n,6} = across_id;
        
    elseif ~isempty(strfind(temp_string,'up'))
        up_index = regexp(temp_string,'\d');
        is_null = regexp(temp_string,'null');
        if ~isempty(is_null)
            up_id = 'null';
        else
            up_id = temp_string(up_index(1):up_index(end));
        end
        base_info{n,7} = up_id;
        
    elseif ~isempty(strfind(temp_string,'down'))      
        down_index = regexp(temp_string,'\d');
        is_null = regexp(temp_string,'null');
        if ~isempty(is_null)
            down_id = 'null';
        else
            down_id = temp_string(down_index(1):down_index(end));
        end
        base_info{n,8} = down_id;
        n = n+1;
    else
    end
       
end

base_info(1,:) = [];

strand_start = strcmp(base_info(:,8),'null'); % find the 3' end
strand_index = find(strand_start==1);  
%%%% note that the base id in the file started from 0

strand_set = {'strand id','base id','postionx','positiony','positionz','basetype','across','up','down'};

for c1 = 1:length(strand_index)
    
    strand_id = sprintf('%d',c1);

    starting_Base_index = strand_index(c1);
       
    base_strand_info = {strand_id,...
                    base_info{starting_Base_index,1},...
                    base_info{starting_Base_index,2},...
                    base_info{starting_Base_index,3},...
                    base_info{starting_Base_index,4},...
                    base_info{starting_Base_index,5},...
                    base_info{starting_Base_index,6},...
                    base_info{starting_Base_index,7},...
                    base_info{starting_Base_index,8}
                    };
    
    strand_set(end+1,:) = base_strand_info;    
    
    down_BaseId = base_info{starting_Base_index,7}; 
    down_Base_index = find(strcmp(base_info(:,1),down_BaseId)==1);
    
    
    while isempty(strfind(down_BaseId,'null')) 
        
        
        down_Base_index = find(strcmp(base_info(:,1),down_BaseId)==1);
        base_strand_info = {strand_id,...
                            base_info{down_Base_index,1},...
                            base_info{down_Base_index,2},...
                            base_info{down_Base_index,3},...
                            base_info{down_Base_index,4},...
                            base_info{down_Base_index,5},...
                            base_info{down_Base_index,6},...
                            base_info{down_Base_index,7},...
                            base_info{down_Base_index,8}
                            };
        strand_set(end+1,:) = base_strand_info;
        
        down_BaseId = base_info{down_Base_index,7};
        
        %fprintf('strand %d: base %s\n',c1,down_BaseId);
                     
    end  
    
end

cd design_file;
writeCSVFilefromStructure('base_strand_info.csv',strand_set);
cd ..;
%%%%% reorder the strand index;
strand_set(1,:) = [];

re_order_strand_set = cell(size(strand_set,1),size(strand_set,2));

new_index = 1:size(strand_set,1);

strand_set_no_strand_index = strand_set;
strand_set_no_strand_index(:,1) = [];

re_order_strand_set_no_strand_index = cell(size(strand_set_no_strand_index,1),size(strand_set_no_strand_index,2));

for c1 = 1:size(re_order_strand_set_no_strand_index,1)
    for c2 = 1:size(re_order_strand_set_no_strand_index,2)
        re_order_strand_set_no_strand_index{c1,c2} = 'null';
        
    end

 end

for c2 = 1:size(strand_set_no_strand_index,1)
    
    original_id = strand_set_no_strand_index{c2,1};
    new_id = num2str(new_index(c2)-1);
    [x,y] = find(strcmp(strand_set_no_strand_index,original_id)==1);
    
    for c3 = 1:length(x)    
        re_order_strand_set_no_strand_index{x(c3),y(c3)} = new_id;
    end
    
end

re_order_strand_set_no_strand_index(:,2) = strand_set_no_strand_index(:,2);
re_order_strand_set_no_strand_index(:,3) = strand_set_no_strand_index(:,3);
re_order_strand_set_no_strand_index(:,4) = strand_set_no_strand_index(:,4);
re_order_strand_set_no_strand_index(:,5) = strand_set_no_strand_index(:,5);

re_order_strand_set(:,1) = strand_set(:,1);
re_order_strand_set(:,2) = re_order_strand_set_no_strand_index(:,1);
re_order_strand_set(:,3) = re_order_strand_set_no_strand_index(:,2);
re_order_strand_set(:,4) = re_order_strand_set_no_strand_index(:,3);
re_order_strand_set(:,5) = re_order_strand_set_no_strand_index(:,4);
re_order_strand_set(:,6) = re_order_strand_set_no_strand_index(:,5);

% cross_base_set =  re_order_strand_set_no_strand_index(:,6);
% 
% for c3 = 1:length(cross_base_set)
%     id = cross_base_set(c3);
%     if isempty(id)
%        cross_base_set(c3) = 'null';
%        fprintf('')
%     else
%     end
%     
% end

re_order_strand_set(:,7) = re_order_strand_set_no_strand_index(:,6);
re_order_strand_set(:,8) = re_order_strand_set_no_strand_index(:,7);
re_order_strand_set(:,9) = re_order_strand_set_no_strand_index(:,8);

strand_set = re_order_strand_set;
scale = 1/0.85;
file_name_conf = sprintf('design_file/init.conf');
fid_conf = fopen(file_name_conf,'w');

file_name_top = sprintf('design_file/sim.top');
fid_top = fopen(file_name_top,'w');

%%%% configuration setup
steps = 1000000000;
max_length = 63;
box_size = max_length*2; % maximum length = max_length
fprintf(fid_conf,'t = %d\n',steps); % define the time steps
fprintf(fid_conf,'b = %d %d %d\n',box_size,box_size,box_size); % define the box sizes
fprintf(fid_conf,'E = 0 0 0\n'); % define the energy

%%%%% topology setup
num_base = size(strand_set,1);
num_strand = length(strand_index);
fprintf(fid_top,'%d ',num_base); % define the box sizes
fprintf(fid_top,'%d\n',num_strand); % define the energy

%for c1 = 26
for c1 = 1:size(strand_set,1)
    strand_id = strand_set{c1,1};
    base_id = strand_set{c1,2};
    posx = str2num(strand_set{c1,3});
    posy = str2num(strand_set{c1,4});
    posz = str2num(strand_set{c1,5});
    
    base_vector = scale * [posx posy posz];
    
    %%% pair base information
    paring_base_id = strand_set{c1,7};
    
    if isempty(strfind(paring_base_id,'null'))
        
       paring_Base_index = find(strcmp(strand_set(:,2),paring_base_id)==1);
       paring_base_posx = str2num(strand_set{paring_Base_index,3});
       paring_base_posy = str2num(strand_set{paring_Base_index,4});
       paring_base_posz = str2num(strand_set{paring_Base_index,5});
       
       paring_base_vector = scale * [paring_base_posx paring_base_posy paring_base_posz];
       
       paring_base5_id = strand_set{paring_Base_index,8}; %% 5'
        if isempty(strfind(paring_base5_id,'null'))

           paring_Base5_index = find(strcmp(strand_set(:,2),paring_base5_id)==1);
           paring_base5_posx = str2num(strand_set{paring_Base5_index,3});
           paring_base5_posy = str2num(strand_set{paring_Base5_index,4});
           paring_base5_posz = str2num(strand_set{paring_Base5_index,5}); 
           
           paring_base5_vector = scale*[paring_base5_posx paring_base5_posy paring_base5_posz];
           
        else
           paring_base5_id = '-1';
           paring_base5_vector = scale*[0 0 1];
           
        end
            
       paring_base3_id = strand_set{paring_Base_index,9}; %% 3'
       
       if isempty(strfind(paring_base3_id,'null'))
        
          paring_Base3_index = find(strcmp(strand_set(:,2),paring_base3_id)==1);
          paring_base3_posx = str2num(strand_set{paring_Base3_index,3});
          paring_base3_posy = str2num(strand_set{paring_Base3_index,4});
          paring_base3_posz = str2num(strand_set{paring_Base3_index,5});          
          paring_base3_vector = scale*[paring_base3_posx paring_base3_posy paring_base3_posz];
       else  
          paring_base3_id = '-1'; 
          paring_base3_vector = [0 0 1];
       end             
    else 
        
        paring_base_id = '-1';
        paring_base_posx = 0;
        paring_base_posy = 0;
        paring_base_posz = 1;
        
        paring_base_vector = [0 0 1];
        
    end
    
    %%% find the vector 2 find the upper base

    up_base_id = strand_set{c1,8};
    if isempty(strfind(up_base_id,'null'))
       up_Base_index = find(strcmp(strand_set(:,2),up_base_id)==1);
       up_base_posx = str2num(strand_set{up_Base_index,3});
       up_base_posy = str2num(strand_set{up_Base_index,4});
       up_base_posz = str2num(strand_set{up_Base_index,5});
       
       up_base_vector = scale*[up_base_posx up_base_posy up_base_posz];
       
       up_base_pairing_id = strand_set{up_Base_index,7};
       
    else
        
        up_base_id = '-1';
        up_base_pairing_id = '-1';
        
        up_base_posx = 0;
        up_base_posy = 0;
        up_base_posz = 1;
        
    end
    
    up_base_vector = scale * [up_base_posx up_base_posy up_base_posz];
   
    down_base_id = strand_set{c1,9};
    if isempty(strfind(down_base_id,'null'))
       down_Base_index = find(strcmp(strand_set(:,2),down_base_id)==1);
       down_base_posx = str2num(strand_set{down_Base_index,3});
       down_base_posy = str2num(strand_set{down_Base_index,4});
       down_base_posz = str2num(strand_set{down_Base_index,5});
       
       down_base_paring_id = strand_set{down_Base_index,7};
       
    else
        %fprintf('true\n')
        down_base_id = '-1';
        down_base_paring_id  = '-1';
        down_base_posx = 0;
        down_base_posy = 0;
        down_base_posz = 1;
    end
    
    down_base_vector = scale*[down_base_posx down_base_posy down_base_posz];
        
    % three backbong vectors, A-base_vector, B-paring_base_vector,
    % A5up_base_vector, A3down_base_vector
    
    backbone_A_to_backbone_B = -base_vector+paring_base_vector;
    backbone_A_to_backbone_B = backbone_A_to_backbone_B/norm(backbone_A_to_backbone_B);
    
    backbone_A_to_posA5_neibor = -base_vector+up_base_vector;
    backbone_A_to_posA5_neibor = backbone_A_to_posA5_neibor/norm(backbone_A_to_posA5_neibor);
    
    backbone_A_to_posA3_neibor = -base_vector+down_base_vector;
    backbone_A_to_posA3_neibor = backbone_A_to_posA3_neibor/norm(backbone_A_to_posA3_neibor);
    
    backboneA_to_backbone_B3 = -base_vector+paring_base3_vector;
    backboneA_to_backbone_B3 = backboneA_to_backbone_B3/norm(backboneA_to_backbone_B3);
    
    backboneA_to_backbone_B5 = -base_vector+paring_base5_vector;
    backboneA_to_backbone_B5 = backboneA_to_backbone_B5/norm(backboneA_to_backbone_B5);
   
    
%     %fprintf('base %s %s; %s %s; %s %s \n',base_id,paring_base_id,...
%                             up_base_id,up_base_pairing_id,...
%                             down_base_id,down_base_paring_id);
    
    
    if str2num(up_base_pairing_id)~=-1 & str2num(paring_base_id)~=-1
        B_5_continuity = abs(str2num(up_base_pairing_id)-str2num(paring_base_id));
    else
        B_5_continuity = 0;
    end
    
    if str2num(down_base_paring_id)~=-1 & str2num(paring_base_id)~=-1
        B_3_continuity = abs(str2num(down_base_paring_id)-str2num(paring_base_id));
    else
        B_3_continuity = 0;
    end

    
    if B_5_continuity==1

       [a1_vector,a3_vector,cm_pos_A] = ...
       Neighbor5_Cal_vector(backbone_A_to_backbone_B,backbone_A_to_posA5_neibor,backboneA_to_backbone_B3);
       cm_pos_A = cm_pos_A+base_vector;
       %fprintf('%d using Neighbor5cal...\n',c1);
       
    elseif B_3_continuity==1
        
       [a1_vector,a3_vector,cm_pos_A] = ...
       Neighbor3_Cal_vector(backbone_A_to_backbone_B,backbone_A_to_posA3_neibor,backboneA_to_backbone_B5);
       cm_pos_A = cm_pos_A+base_vector;
       %fprintf('%d using Neighbor3cal...\n',c1);
       
    else
        
       %fprintf('%d using not using...!!!!!!!\n',c1);
       cm_pos_A = base_vector;
       a3_vector = a3_vector;
       a1_vector = a1_vector;

    end
    
    base = strand_set{c1,6};
    base = base(1);
    
    
    %fprintf('%d %d:\n',c1,a1_vector(1)^2+a1_vector(2)^2+a1_vector(3)^2);
    
    fprintf(fid_top,'%s %s %s %s \n',strand_id,base,down_base_id,up_base_id); % position
    
    fprintf(fid_conf,'%d %d %d ',cm_pos_A); % position
    fprintf(fid_conf,'%d %d %d ',a1_vector); % vector1
    fprintf(fid_conf,'%d %d %d ',a3_vector); % vector2
    fprintf(fid_conf,'0 0 0 '); % velocity
    fprintf(fid_conf,'0 0 0\n'); % angle momentum

      
end

%open sim.top

%system('~/old_oxdna/oxDNA/UTILS/traj2chimera.py init.conf sim.top');
%system('/Applications/Chimera.app/Contents/MacOS/chimera init.conf.pdb chimera.com');
% system('cp ~/old_oxdna/input .');
% system('cp ~/old_oxdna/input_semirelax .');
% system('~/old_oxdna/oxDNA/bin/oxDNA input');
% system('~/old_oxdna/oxDNA/bin/oxDNA input_semirelax');












