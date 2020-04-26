package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.getHighestMatchType;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionWriter
{
    private final IsofoxConfig mConfig;
    private BufferedWriter mReadWriter;
    private BufferedWriter mFragmentWriter;

    public FusionWriter(final IsofoxConfig config)
    {
        mConfig = config;

        mReadWriter = null;
        mFragmentWriter = null;
    }

    public void close()
    {
        closeBufferedWriter(mReadWriter);
        closeBufferedWriter(mFragmentWriter);
    }

    public static final String FUSION_FILE_ID = "fusions.csv";

    public void writeFusionData(Map<String,List<FusionReadData>> fusionCandidates)
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile(FUSION_FILE_ID);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write(FusionReadData.csvHeader());
            writer.newLine();

            for(Map.Entry<String, List<FusionReadData>> entry : fusionCandidates.entrySet())
            {
                for (final FusionReadData fusion : entry.getValue())
                {
                    writer.write(fusion.toCsv());
                    writer.newLine();
                }
            }

            writer.close();

            if(mConfig.WriteChimericReads)
            {
                for (List<FusionReadData> fusions : fusionCandidates.values())
                {
                    for (FusionReadData fusion : fusions)
                    {
                        for (FusionFragment fragment : fusion.getAllFragments())
                        {
                            writeFragmentData(fragment, fusionId(fusion.id()));
                            writeReadData(fragment.getReads(), fusionId(fusion.id()));
                        }
                    }
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fusions file: {}", e.toString());
        }
    }

    public void writeReadData(final List<ReadRecord> reads, final String groupStatus)
    {
        if(!mConfig.WriteChimericReads)
            return;

        try
        {
            if(mReadWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("chimeric_reads.csv");

                mReadWriter = createBufferedWriter(outputFileName, false);
                mReadWriter.write("ReadSetCount,ReadId,FusionGroup,Chromosome,PosStart,PosEnd,Cigar,InsertSize");
                mReadWriter.write(",Secondary,Supplementary,NegStrand,ProperPair,SuppAlign,TransExons,BestMatch,TransExonData");
                mReadWriter.newLine();
            }

            for(final ReadRecord read : reads)
            {
                mReadWriter.write(String.format("%s,%s,%s,%s,%d,%d,%s,%d",
                        reads.size(), read.Id, groupStatus, read.Chromosome,
                        read.PosStart, read.PosEnd, read.Cigar.toString(), read.fragmentInsertSize()));

                mReadWriter.write(String.format(",%s,%s,%s,%s,%s",
                        read.isSecondaryAlignment(), read.isSupplementaryAlignment(), read.isNegStrand(), read.isProperPair(),
                        read.getSuppAlignment() != null ? read.getSuppAlignment().replaceAll(",", ";") : "NONE"));

                // log the transcript exons affected, and the highest matching transcript
                String transExonData = "";

                int transExonRefCount = 0;
                RegionMatchType highestTransMatchType = getHighestMatchType(read.getTransExonRefs().keySet());

                if(highestTransMatchType != NONE)
                {
                    for (final TransExonRef transExonRef : read.getTransExonRefs().get(highestTransMatchType))
                    {
                        ++transExonRefCount;
                        transExonData = appendStr(transExonData, String.format("%s:%d",
                                transExonRef.TransName, transExonRef.ExonRank), ';');
                    }
                }

                mReadWriter.write(String.format(",%d,%s,%s",
                        transExonRefCount, highestTransMatchType, transExonRefCount == 0 ? "NONE" : transExonData));
                mReadWriter.newLine();
            }

        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return;
        }

    }

    public void writeUnfusedFragments(final Map<String,List<FusionFragment>> unfusedFragments)
    {
        if(!mConfig.WriteChimericReads)
            return;

        for(List<FusionFragment> fragments : unfusedFragments.values())
        {
            for(FusionFragment fragment : fragments)
            {
                writeFragmentData(fragment, "UNFUSED");
                writeReadData(fragment.getReads(), "UNFUSED");
            }
        }
    }

    private void writeFragmentData(final FusionFragment fragment, final String fusionId)
    {
        if(!mConfig.WriteChimericReads)
            return;

        try
        {
            if(mFragmentWriter == null)
            {
                final String outputFileName = mConfig.formOutputFile("chimeric_frags.csv");

                mFragmentWriter = createBufferedWriter(outputFileName, false);
                mFragmentWriter.write("ReadId,ReadCount,FusionGroup,Type,SameGene,ScCount");

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    final String prefix = se == SE_START ? "Start" : "End";
                    mFragmentWriter.write(",Chr" + prefix);
                    mFragmentWriter.write(",JuncPos" + prefix);
                    mFragmentWriter.write(",JuncOrient" + prefix);
                    mFragmentWriter.write(",Region" + prefix);
                    mFragmentWriter.write(",JuncType" + prefix);
                }

                mFragmentWriter.newLine();
            }

            mFragmentWriter.write(String.format("%s,%d,%s,%s,%s,%d",
                    fragment.readId(), fragment.getReads().size(), fusionId, fragment.type(),
                    fragment.geneCollections()[SE_START] == fragment.geneCollections()[SE_END],
                    fragment.getReads().stream().filter(x -> x.containsSoftClipping()).count()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                mFragmentWriter.write(String.format(",%s,%d,%d,%s,%s",
                        fragment.chromosomes()[se], fragment.junctionPositions()[se], fragment.junctionOrientations()[se],
                        fragment.regionMatchTypes()[se], fragment.junctionTypes()[se]));
            }

            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return;
        }

    }

}
