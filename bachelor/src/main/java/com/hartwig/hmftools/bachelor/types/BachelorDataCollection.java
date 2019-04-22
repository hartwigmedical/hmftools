package com.hartwig.hmftools.bachelor.types;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class BachelorDataCollection
{
    // SampleId,Program,Id,Gene,TranscriptId,Chromosome,Position,Ref,Alt
    // CodingEffect,Effect,Annotations,HgvsProtein,IsHomozygous,PhredScore,HgvsCoding
    // MatchType,HasDepthInfo,GermlineAltCount,GermlineReadDepth,TumorAltCount,TumorReadDepth,CodonInfo
    private static final int COL_INDEX_SAMPLE = 0;
    private static final int COL_INDEX_PROGRAM = 1;
    private static final int COL_INDEX_SV_ID = 2;
    private static final int COL_INDEX_GENE = 3;
    private static final int COL_INDEX_TRAN_ID = 4;
    private static final int COL_INDEX_CHR = 5;
    private static final int COL_INDEX_POS = 6;
    private static final int COL_INDEX_REF = 7;
    private static final int COL_INDEX_ALT = 8;

    private static final int COL_INDEX_CODING_EFFECT = 9;
    private static final int COL_INDEX_EFFECTS = 10;
    private static final int COL_INDEX_ANNOTS = 11;
    private static final int COL_INDEX_PROTEIN = 12;
    private static final int COL_INDEX_HZ = 13;
    private static final int COL_INDEX_PHRED = 14;
    private static final int COL_INDEX_CODING = 15;
    private static final int COL_INDEX_MATCH_TYPE = 16;
    private static final int COL_INDEX_HAS_DEPTH = 17;
    private static final int COL_INDEX_GL_ALT_COUNT = 18;
    private static final int COL_INDEX_GL_READ_DEPTH = 19;
    private static final int COL_INDEX_TUMOR_ALT_COUNT = 20;
    private static final int COL_INDEX_TUMOR_READ_DEPTH = 21;
    private static final int COL_INDEX_CODON_INFO = 22;
    private static final int COL_INDEX_CLINVAR_MATCH = 23;
    private static final int COL_INDEX_CLINVAR_SIG = 24;
    private static final int COL_INDEX_CLINVAR_SIG_INFO = 25;

    private static final int BACHELOR_CSV_FIELD_COUNT = COL_INDEX_CODON_INFO + 1;

    private static final Logger LOGGER = LogManager.getLogger(BachelorDataCollection.class);

    private List<String> mLimitedSampleList;
    private List<BachelorGermlineVariant> mGermlineVariants;
    private int mFileIndex;
    private int mMaxReadCount;

    private BufferedReader mFileReader;

    public BachelorDataCollection()
    {
        mGermlineVariants = Lists.newArrayList();
        mFileIndex = 0;
        mMaxReadCount = 0;
        mFileReader = null;
    }

    public void setSampleIds(final List<String> sampleIdsList)
    {
        mLimitedSampleList.addAll(sampleIdsList);
    }
    public void setMaxReadCount(int maxReadCount) { mMaxReadCount = maxReadCount; }

    public boolean loadBachelorData(final String filename)
    {
        if (filename.isEmpty())
            return false;

        try
        {
            mFileReader = new BufferedReader(new FileReader(filename));

            mFileReader.readLine(); // skip header
        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read bachelor input CSV file({})", filename);
            return false;
        }

        return true;
    }

    public boolean processBachelorData()
    {
        mGermlineVariants.clear();

        try
        {
            String line = null;
            boolean exitNextSample = false;
            String currentSampleId = "";

            while ((line = mFileReader.readLine()) != null)
            {
                if(line.isEmpty())
                    break;

                ++mFileIndex;

                // parse CSV data
                String[] items = line.split(",", -1);

                if (items.length < BACHELOR_CSV_FIELD_COUNT)
                {
                    LOGGER.warn("sample({}) invalid item count({}), fileIndex({})",
                            currentSampleId, items.length, mFileIndex);
                    continue;
                }

                final String sampleId = items[COL_INDEX_SAMPLE];

                if(!mLimitedSampleList.isEmpty() && !mLimitedSampleList.contains(sampleId))
                    continue;

                if(!sampleId.equals(currentSampleId))
                {
                    currentSampleId = sampleId;

                    if(exitNextSample)
                    {
                        LOGGER.info("halting read at {} records, fileIndex({})", mGermlineVariants.size(), mFileIndex);
                        return true;
                    }
                }

                // check for annotations with ',' which impacts string splitting
                if(items.length > BACHELOR_CSV_FIELD_COUNT)
                {
                    checkAnnotationItems(items);
                }

                try
                {

                    BachelorGermlineVariant bachRecord = new BachelorGermlineVariant(sampleId,
                            items[COL_INDEX_PROGRAM],
                            items[COL_INDEX_SV_ID],
                            items[COL_INDEX_GENE],
                            items[COL_INDEX_TRAN_ID],
                            items[COL_INDEX_CHR],
                            Long.parseLong(items[COL_INDEX_POS]),
                            items[COL_INDEX_REF],
                            items[COL_INDEX_ALT],
                            CodingEffect.valueOf(items[COL_INDEX_CODING_EFFECT]),
                            items[COL_INDEX_EFFECTS],
                            items[COL_INDEX_ANNOTS],
                            items[COL_INDEX_PROTEIN],
                            Boolean.parseBoolean(items[COL_INDEX_HZ]),
                            Integer.parseInt(items[COL_INDEX_PHRED]),
                            items[COL_INDEX_CODING],
                            items[COL_INDEX_MATCH_TYPE],
                            items[COL_INDEX_CODON_INFO],
                            Boolean.parseBoolean(items[COL_INDEX_CLINVAR_MATCH]),
                            items[COL_INDEX_CLINVAR_SIG],
                            items[COL_INDEX_CLINVAR_SIG_INFO]);

                    boolean hasDepthInfo = Boolean.parseBoolean(items[COL_INDEX_HAS_DEPTH]);

                    int glAltCount = Integer.parseInt(items[COL_INDEX_GL_ALT_COUNT]);
                    int glReadDepth = Integer.parseInt(items[COL_INDEX_GL_READ_DEPTH]);
                    bachRecord.setGermlineData(glAltCount, glReadDepth);

                    if (hasDepthInfo)
                    {
                        int tumorAltCount = Integer.parseInt(items[COL_INDEX_TUMOR_ALT_COUNT]);
                        int tumorReadDepth = Integer.parseInt(items[COL_INDEX_TUMOR_READ_DEPTH]);

                        bachRecord.setReadData(tumorAltCount, tumorReadDepth);
                    }

                    mGermlineVariants.add(bachRecord);
                }
                catch(Exception nfe)
                {
                    LOGGER.debug("line parse error({}) fileIndex({}) line: {}", nfe.toString(), mFileIndex, line);
                    return false;
                }

                if(!exitNextSample && mMaxReadCount > 0 && mGermlineVariants.size() >= mMaxReadCount)
                    exitNextSample = true;
            }

            LOGGER.debug("loaded {} bachelor records", mGermlineVariants.size());
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read bachelor input CSV file: {}", e.toString());
            return false;
        }

        return !mGermlineVariants.isEmpty();
    }

    private static void checkAnnotationItems(String[] items)
    {
        String extraAnnots = "";
        boolean hasExtraAnnots = false;
        int extraItems = 0;

        for (int i = COL_INDEX_ANNOTS + 1; i < items.length; ++i)
        {
            extraAnnots += "," + items[i];
            ++extraItems;

            if (items[i].contains("||"))
            {
                hasExtraAnnots = true;
                break;
            }
        }

        if (!hasExtraAnnots)
            return;

        // otherwise shift them back down
        items[COL_INDEX_ANNOTS] += extraAnnots;

        for (int i = COL_INDEX_ANNOTS + 1; i <= COL_INDEX_CODON_INFO; ++i)
        {
            items[i] = items[i + extraItems];
        }
    }

    public final List<BachelorGermlineVariant> getBachelorVariants() { return mGermlineVariants; }

}

