package com.hartwig.hmftools.cup.sigs;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.sigs.RefSignatures.populateRefSigContributions;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.DataUtils;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.utils.GenericDataCollection;
import com.hartwig.hmftools.common.utils.GenericDataLoader;
import com.hartwig.hmftools.cup.svs.SvData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class SignatureDataLoader
{
    public static SigMatrix loadSampleCountsFromCohortFile(final String filename, final Map<String,Integer> sampleCountsIndex)
    {
        final GenericDataCollection collection = GenericDataLoader.loadFile(filename);

        for(int s = 0; s < collection.getFieldNames().size(); ++s)
        {
            final String sampleId = collection.getFieldNames().get(s);
            sampleCountsIndex.put(sampleId, s);
        }

        SigMatrix sampleCounts = DataUtils.createMatrixFromListData(collection.getData());
        sampleCounts.cacheTranspose();

        return sampleCounts;
    }


    public static void loadSampleCountsFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String, SvData> sampleSvData)
    {
        if(dbAccess == null)
            return;

        for(final String sampleId : sampleIds)
        {
            final List<SignatureAllocation> sigAllocations = dbAccess.readSignatureAllocations(sampleId);

        }
    }

    public static void loadSigContribsFromCohortFile(final String filename, final Map<String,Map<String,Double>> sampleSigContributions)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                // SampleId,SigName,SigContrib,SigPercent
                final String[] items = line.split(DATA_DELIM, -1);
                String sampleId = items[0];
                String sigName = items[1];
                double sigContrib = Double.parseDouble(items[2]);

                Map<String,Double> sigContribs = sampleSigContributions.get(sampleId);

                if(sigContribs == null)
                {
                    sigContribs = Maps.newHashMap();
                    sampleSigContributions.put(sampleId, sigContribs);
                }

                sigContribs.put(sigName, sigContrib);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sig contribution data file({}): {}", filename, e.toString());
        }
    }

    public static void loadRefSigContribPercentiles(final String filename, final Map<String,Map<String,double[]>> refCancerSigContribPercentiles)
    {
        if(filename.isEmpty())
            return;

        populateRefSigContributions(filename, refCancerSigContribPercentiles);
    }

    public static boolean loadRefSampleCounts(final String filename, SigMatrix refSampleCounts, final List<String> refSampleNames)
    {
        final GenericDataCollection collection = GenericDataLoader.loadFile(filename);

        refSampleNames.addAll(collection.getFieldNames());
        refSampleCounts = DataUtils.createMatrixFromListData(collection.getData());
        refSampleCounts.cacheTranspose();

        return true;
    }


}
