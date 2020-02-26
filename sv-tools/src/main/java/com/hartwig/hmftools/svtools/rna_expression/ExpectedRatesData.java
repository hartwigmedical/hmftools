package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.RE_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sig_analyser.common.SigMatrix;

public class ExpectedRatesData
{
    public final String GeneId;

    // equivalent of buckets - 0-N transcripts and the fragment type (eg SHORT, SPLICED etc)
    public final List<String> Categories;

    // equivalent of signature names - all transcript names (ie StableId) and an UNSPLICED definition
    public final List<String> TranscriptIds;

    private SigMatrix mTranscriptDefinitions;

    public ExpectedRatesData(final String geneId)
    {
        GeneId = geneId;
        Categories = Lists.newArrayList();
        TranscriptIds = Lists.newArrayList();
        mTranscriptDefinitions = null;
    }

    public SigMatrix getTranscriptDefinitions() { return mTranscriptDefinitions; }


    public boolean validData()
    {
        if(Categories.isEmpty() || mTranscriptDefinitions == null)
            return false;

        if(mTranscriptDefinitions.Cols != TranscriptIds.size())
            return false;

        if(mTranscriptDefinitions.Rows != Categories.size())
            return false;

        return true;
    }

    public void initialiseTranscriptDefinitions()
    {
        if(Categories.isEmpty() || TranscriptIds.isEmpty())
            return;

        mTranscriptDefinitions = new SigMatrix(Categories.size(), TranscriptIds.size());
    }

    // GeneId,GeneName,TransId,Category,Rate
    public static final int ER_COL_GENE_ID = 0;
    public static final int ER_COL_GENE_NAME = 1;
    public static final int ER_COL_TRANS_NAME = 2;
    public static final int ER_COL_CAT = 3;
    public static final int ER_COL_RATE = 4;

    public void buildDefinitionsFromFileData(final List<String[]> geneRatesData)
    {
        initialiseTranscriptDefinitions();

        double[][] matrixData = mTranscriptDefinitions.getData();

        for(String[] data : geneRatesData)
        {
            String transName = data[ER_COL_TRANS_NAME];
            String categoryStr = data[ER_COL_CAT];

            int transIndex = getTranscriptIndex(transName);
            int catIndex = getCategoryIndex(categoryStr);
            double rate = Double.parseDouble(data[ER_COL_RATE]);
            matrixData[catIndex][transIndex] = rate;
        }
    }

    private int getTranscriptIndex(final String trans)
    {
        for(int i = 0; i < TranscriptIds.size(); ++i)
        {
            if(TranscriptIds.get(i).equals(trans))
                return i;
        }

        return -1;
    }

    public int getCategoryIndex(final String category)
    {
        for(int i = 0; i < Categories.size(); ++i)
        {
            if(Categories.get(i).equals(category))
                return i;
        }

        return -1;
    }

    public void addCategory(final String category)
    {
        if(!Categories.contains(category))
            Categories.add(category);
    }

}
