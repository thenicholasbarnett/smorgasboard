#include <iostream>
#include <iostream>
void Newtxtfile(){

// naming output and input txt file
    int i = 0;
    string inputtxtfile = Form("2024ppRef_MC_%d_filenames.txt",i);
    string outputtxtfile = Form("2024ppRef_MC_%d_filenames_CRAB.txt",i);

// opening input txt file
    ifstream myfile(inputtxtfile);

// making a string for each line in input txt file
    string filename;

// creating output txt file
    std::ofstream outputFile(outputtxtfile);

// while looping over each line in input text file
    while(getline(myfile, filename)){

// taking out first nine characters in each string composed of each line in input txt file
        filename.erase(0,8);

// writing string that is the remainder of each line in input txt file into the output txt file
        outputFile << filename+"\n";
    }

// closing the output txt file
    outputFile.close();
}
