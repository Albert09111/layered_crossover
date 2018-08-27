function position = get_position(index)
    id = str2num(index);
    base_info = loadCSVabsFileRawtoStructure('design_file/base_strand_info.csv','');
    base_x = str2num(base_info{id+1,3});
    base_y = str2num(base_info{id+1,4});
    base_z = str2num(base_info{id+1,5});
    position = [base_x,base_y,base_z];
end