package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class IndelAnnotator
{
    private final DatabaseAccess mDbAccess;

    private final List<IndelData> mSampleIndelData;

    private static final Logger LOGGER = LogManager.getLogger(IndelAnnotator.class);

    public IndelAnnotator(final DatabaseAccess dbAccess)
    {
        mDbAccess = dbAccess;
        mSampleIndelData = Lists.newArrayList();
    }
    public final List<IndelData> getSampleIndelData() { return mSampleIndelData; }

    public void loadIndels(final String sampleId)
    {
        LOGGER.info("retrieving SNV data for sample({})", sampleId);

        final List<SomaticVariant> variants = mDbAccess.readSomaticVariants(sampleId, VariantType.INDEL);

        LOGGER.info("sample({}) processing {} variants", sampleId, variants.size());

        for(final SomaticVariant variant : variants)
        {
            if (variant.isFiltered())
                continue;

            IndelData indel = new IndelData(variant.chromosome(), variant.position(), variant.ref(), variant.alt(),
                    variant.microhomology(), variant.repeatSequence(), variant.repeatCount(), variant.ploidy());

            // which ones to keep or only decide after analysing profile of the sample?
        }
    }

}
