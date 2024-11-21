#include <iostream>
#include <iostream>
void Newtxtfile(){
    string oldtxtfile = Form("2024ppRef_MC_%d_filenames.txt",i);
    string newtxtfile = Form("2024ppRef_MC_%d_filenames_CRAB.txt",i);
    ifstream myfile(oldtxtfile);
    string filename;
    std::ofstream outputFile(newtxtfile);
    if (outputFile.is_open()){
        while(getline(myfile, filename)){
            filename.erase(0,8);
            outputFile << filename+"\n";
        }
        outputFile.close();
    }
}
