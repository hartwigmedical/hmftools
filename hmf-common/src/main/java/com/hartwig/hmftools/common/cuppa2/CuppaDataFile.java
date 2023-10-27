package com.hartwig.hmftools.common.cuppa2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

public class CuppaDataFile
{
    public final String Filename;
    public final List<CuppaPrediction> CuppaPredictions;
    public final boolean HasRnaData;
    public final String MainCombinedClassifier;

    CuppaDataFile(final String filename) throws IOException {
        Filename = filename;
        CuppaPredictions = readTable(filename);
        HasRnaData = checkHasRnaData();
        MainCombinedClassifier = getMainCombinedClassifier();
    }

    // Read table --------------------------------
    private static double parseDouble(String string) {
        if(string.length()==0) {
            string = "NaN";
        } else if (string.equals("inf")) {
            string = "Infinity";
        } else if (string.equals("-inf")) {
            string = "-Infinity";
        }

        return Double.parseDouble(string);
    }

    private static String parseString(final String string){
        if(string.length() > 0) {
            return string;
        }
        return "none";
    }

    private static List<CuppaPrediction> readTable(final String filename) throws IOException
    {
        String delimiter = "\t";

        BufferedReader fileReader = new BufferedReader(new FileReader(filename));

        String line = fileReader.readLine();
        final Map<String, Integer> fieldsMap = createFieldsIndexMap(line, delimiter);
        int sampleIdIndex = fieldsMap.get(CuppaPrediction.FLD_SAMPLE_ID);
        int dataTypeIndex = fieldsMap.get(CuppaPrediction.FLD_DATA_TYPE);
        int clfGroupIndex = fieldsMap.get(CuppaPrediction.FLD_CLF_GROUP);
        int clfNameIndex = fieldsMap.get(CuppaPrediction.FLD_CLF_NAME);
        int featNameIndex = fieldsMap.get(CuppaPrediction.FLD_FEAT_NAME);
        int featValueIndex = fieldsMap.get(CuppaPrediction.FLD_FEAT_VALUE);
        int cancerTypeIndex = fieldsMap.get(CuppaPrediction.FLD_CANCER_TYPE);
        int dataValueIndex = fieldsMap.get(CuppaPrediction.FLD_DATA_VALUE);

        List<CuppaPrediction> cuppaPredictions = new ArrayList<>();
        while((line = fileReader.readLine()) != null)
        {
            String[] rowValues = line.split(delimiter, -1);

            String sampleId = parseString(rowValues[sampleIdIndex]);

            String dataTypeStr = parseString(rowValues[dataTypeIndex]);
            Categories.DataType dataType = Categories.DataType.valueOf(dataTypeStr);

            String clfGroupStr = parseString(rowValues[clfGroupIndex]);
            Categories.ClfGroup clfGroup = Categories.ClfGroup.valueOf(clfGroupStr);

            String clfNameStr = parseString(rowValues[clfNameIndex]);
            Categories.ClfName clfName = Categories.ClfName.valueOf(clfNameStr);

            String featName = parseString(rowValues[featNameIndex]);
            double featValue = parseDouble(rowValues[featValueIndex]);
            String cancerType = parseString(rowValues[cancerTypeIndex]);
            double dataValue = parseDouble(rowValues[dataValueIndex]);

            CuppaPrediction cuppaPrediction = new CuppaPrediction(
                    sampleId, dataType, clfGroup, clfName,
                    featName, featValue, cancerType, dataValue
            );

            cuppaPredictions.add(cuppaPrediction);
        }

        return cuppaPredictions;
    }

    // --------------------------------
    private boolean checkHasRnaData(){
        for(CuppaPrediction cuppaPrediction : CuppaPredictions){
            if(!cuppaPrediction.DataType.equals(Categories.DataType.prob)){
                continue;
            }

            if(cuppaPrediction.ClfName.equals(Categories.ClfName.rna_combined) & !Double.isNaN(cuppaPrediction.DataValue)){
                return true;
            }
        }

        return false;
    }

    private String getMainCombinedClassifier(){
        if(HasRnaData){ return "combined"; }
        return "dna_combined";
    }

    // --------------------------------


    // Misc --------------------------------
    public void printPredictions(int n){
        int i = 0;
        while(i < n){
            System.out.println(CuppaPredictions.get(i).toString());
            i++;
        }
    }

    public static void main(String[] args) throws IOException {
        String filename = "/Users/lnguyen/Desktop/cuppa_vis_data.tsv";
        CuppaDataFile cuppaDataFile = new CuppaDataFile(filename);



        // cuppaDataFile.printPredictions(10);
        // System.out.println(cuppaDataFile.mHasRnaData);
        //System.out.println(cuppaDataFile.CuppaPredictions.get(0).DataType);

    }
}


