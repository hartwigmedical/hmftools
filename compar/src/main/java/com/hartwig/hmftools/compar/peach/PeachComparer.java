package com.hartwig.hmftools.compar.peach;

import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.Category.PEACH;
import static com.hartwig.hmftools.compar.common.CommonUtils.fileExists;
import static com.hartwig.hmftools.compar.peach.PeachData.FLD_ALLELE_COUNT;
import static com.hartwig.hmftools.compar.peach.PeachData.FLD_DRUGS;
import static com.hartwig.hmftools.compar.peach.PeachData.FLD_FUNCTION;
import static com.hartwig.hmftools.compar.peach.PeachData.FLD_PRESCRIPTION_URLS;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.SampleFileSources;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class PeachComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public PeachComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category()
    {
        return PEACH;
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_ALLELE_COUNT, FLD_FUNCTION, FLD_DRUGS, FLD_PRESCRIPTION_URLS);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        List<PeachGenotype> genotypes = dbAccess.readPeachGenotypes(sampleId);
        return genotypes.stream().map(g -> new PeachData(g)).collect(Collectors.toList());
    }

    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final SampleFileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        String fileName = determineFileName(sampleId, germlineSampleId, fileSources);
        try
        {
            PeachGenotypeFile.read(fileName).forEach(g -> comparableItems.add(new PeachData(g)));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load Peach data: {}", sampleId, e.toString());
            return null;
        }
        return comparableItems;
    }

    private static String determineFileName(final String sampleId, final String germlineSampleId, final SampleFileSources fileSources)
    {
        final String currentFileName = PeachGenotypeFile.generateFileName(fileSources.peach(), germlineSampleId);
        final String oldFileName = PeachGenotypeFile.generateOldPythonFileName(fileSources.peach(), sampleId);
        if(!fileExists(currentFileName) && fileExists(oldFileName))
        {
            return oldFileName;
        }
        else
        {
            return currentFileName;
        }
    }
}
