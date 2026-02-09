package com.hartwig.hmftools.isofox.expression;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.expression.cohort.GeneratePanelNormalisation.FLD_TPM_ADJUST_FACTOR;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.results.GeneResult;

public class PanelTpmNormaliser
{
    private final Map<String,Double> mGeneNormFactors;
    private final boolean mEnabled;

    public PanelTpmNormaliser(final String filename)
    {
        mGeneNormFactors = Maps.newHashMap();
        loadNormalisations(filename);
        mEnabled = !mGeneNormFactors.isEmpty();
    }

    public void applyNormalisation(final List<GeneCollectionSummary> allGeneSummaries)
    {
        if(!mEnabled)
            return;

        for(final GeneCollectionSummary geneSummary : allGeneSummaries)
        {
            for(GeneResult geneResult : geneSummary.GeneResults)
            {
                Double adjustFactor = mGeneNormFactors.get(geneResult.Gene.GeneId);

                if(adjustFactor == null || adjustFactor <= 0)
                {
                    ISF_LOGGER.warn("gene({}:{}) missing or invalid panel TPM adjustment factor({})",
                            geneResult.Gene.GeneId, geneResult.Gene.GeneName, adjustFactor);
                }
                else
                {
                    geneResult.applyTpmAdjustFactor(adjustFactor);
                }
            }
        }
    }

    private void loadNormalisations(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String fileDelim = inferFileDelimiter(filename);

            Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), fileDelim);
            lines.remove(0);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int adjustIndex = fieldsMap.get(FLD_TPM_ADJUST_FACTOR);

            for(final String data : lines)
            {
                String[] values = data.split(fileDelim);

                String geneId = values[geneIdIndex];
                double adjustFactor = Double.parseDouble(values[adjustIndex]);
                mGeneNormFactors.put(geneId, adjustFactor);

            }

            ISF_LOGGER.info("loaded {} panel gene TPM adjustment factors from file({})", mGeneNormFactors.size(), filename);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load panel gene TPM adjustment factors from file({}): {}", filename, e.toString());
        }
    }
}
