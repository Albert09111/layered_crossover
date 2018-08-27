function array = loadCSVabsFileRawtoStructure(FileName, PathName)
%function array = loadCSVabsFileRawtoStructure(FileName, PathName)

%[FileName,PathName] = uigetfile('*.csv','Select an .csv file corresponding to the Wavelength; Absorbance matrix...');
TotalName = [PathName,FileName];
FileID = fopen(TotalName);

%load all lines of file
array = [];
counter = 1;
temp_string = fgetl(FileID);
while length(temp_string) ~= 0
    if isnumeric(temp_string) == 1
        break;
    end
    line1 = parse_csv_text_string(temp_string);
    for c1 = 1:length(line1)
        array{counter,c1} = line1{c1};
    end
    temp_string = fgetl(FileID);
    counter = counter + 1;
end
fclose(FileID);