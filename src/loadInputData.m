userInput = fopen(inFile,'r');
while ~feof(userInput)
    eval(fgetl(userInput))
end
fclose(userInput);