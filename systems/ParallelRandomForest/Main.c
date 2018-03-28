
#include "ArffImporter.h"
#include "Classifier.h"

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

int main( int argc, char *argv[])
{
    if (argc < 3) {
	printf("You need to give three parameters: training file, testing file, andoutput file\n");
	return 1;
    }
    //get the train filename

    char* trainFileName;
    char * argFilename1 = getCmdOption(argv, argv + argc, "-trainfile");
    if(argFilename1)
    {
        trainFileName = argFilename1;
	printf("Training file is: %s\n", trainFileName);
    }

    char* testFileName;
    char * argFilename2 = getCmdOption(argv, argv + argc, "-testfile");
    if(argFilename2)
    {
        testFileName = argFilename2;
	printf("Testing file is: %s\n", testFileName);
    }

    char* resultFileName;
    char * argFilename3 = getCmdOption(argv, argv + argc, "-resultfile");
    if(argFilename3)
    {
        resultFileName = argFilename3;
	printf("Output file is: %s\n", resultFileName);
    }

    ArffImporter trainSetImporter;
    trainSetImporter.Read(trainFileName);

    ArffImporter testSetImporter;
    testSetImporter.Read(testFileName);

    Classifier classifier;
    classifier.Train(
        trainSetImporter.GetInstances(), 
        trainSetImporter.GetFeatures(), 
        trainSetImporter.GetClassAttr(),
        trainSetImporter.GetNumInstances() );
    classifier.Classify(
        testSetImporter.GetInstances(),
        testSetImporter.GetNumInstances(),
        resultFileName );

    return 0;
}
