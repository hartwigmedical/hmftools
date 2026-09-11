package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.rna.RnaFusionFile.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.NONE;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.RegionMatchType;

public class FusionWriter
{
    private final IsofoxConfig mConfig;
    private BufferedWriter mFusionWriter;
    private BufferedWriter mFusionReadWriter;
    private BufferedWriter mPassingFusionWriter;
    private BufferedWriter mFragmentWriter;
    private final boolean mWriteReads;
    private final boolean mWriteFragments;

    private int mNextFusionId;

    public static final String UNFILTERED_FUSION_FILE_ID = "fusions.tsv";

    public FusionWriter(final IsofoxConfig config)
    {
        mConfig = config;
        mWriteReads = mConfig.Fusions.WriteChimericReads;
        mWriteFragments = mConfig.Fusions.WriteChimericFragments;

        mFusionWriter = null;
        mFragmentWriter = null;
        mFusionReadWriter = null;
        mNextFusionId = 0;

        initialiseFusionWriters();
        initialiseFragmentWriter();
    }

    public synchronized int getNextFusionId() { return mNextFusionId++; }

    public void close()
    {
        closeBufferedWriter(mFusionWriter);
        closeBufferedWriter(mPassingFusionWriter);
        closeBufferedWriter(mFragmentWriter);
        closeBufferedWriter(mFusionReadWriter);
    }

    private void initialiseFusionWriters()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            mFusionWriter = createBufferedWriter(mConfig.formOutputFile(UNFILTERED_FUSION_FILE_ID), false);
            mFusionWriter.write(FusionData.header());
            mFusionWriter.newLine();

            mPassingFusionWriter = createBufferedWriter(mConfig.formOutputFile(PASS_FUSION_FILE_ID), false);
            mPassingFusionWriter.write(RnaFusionFile.header());
            mPassingFusionWriter.newLine();

            mFusionReadWriter = initialiseReadWriter();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to create fusions file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseReadWriter()
    {
        // TODO: purge if too similar to chimeric read writing, or update to TSV format
        try
        {
            final String outputFileName = mConfig.formOutputFile("fusion_reads.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("ReadGroupCount,ReadId,FusionGroup,Chromosome,PosStart,PosEnd,Orientation,Cigar");
            writer.write(",Flags,HasSupplAlign,SuppData,BasesStart,BasesEnd,MateChr,MatePosStart");
            writer.write(",GeneSetStart,GeneSetEnd,GenicStart,GenicEnd,InterGeneSplit");
            writer.write(",MappedCoords,ScRegionsMatchedStart,ScRegionsMatchedEnd");
            writer.write(",TopTransMatch,TransExonData,UpperTopTransMatch,UpperTransExonData");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeFusionData(
            final List<FusionData> fusions, final List<FusionData> passingFusions, final Map<String,List<FusionReadData>> fusionCandidates)
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            for(FusionData fusionData : fusions)
            {
                mFusionWriter.write(fusionData.toTsv());
                mFusionWriter.newLine();
            }

            for(FusionData fusionData : passingFusions)
            {
                RnaFusion fusion = fusionData.buildRnaFusion();
                mPassingFusionWriter.write(RnaFusionFile.write(fusion));
                mPassingFusionWriter.newLine();
            }

            if(mWriteReads || mWriteFragments)
            {
                for(List<FusionReadData> fusionCandidate : fusionCandidates.values())
                {
                    for(FusionReadData fusion : fusionCandidate)
                    {
                        for(List<FusionFragment> fragments : fusion.getFragments().values())
                        {
                            for(FusionFragment fragment : fragments)
                            {
                                if(mWriteFragments)
                                    writeFragmentData(fragment, fusionId(fusion.id()));

                                if(mWriteReads)
                                    writeReadData(fragment.readId(), fragment.reads(), fusionId(fusion.id()));
                            }
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

    public synchronized void writeUnfusedFragments(final List<FusionFragment> fragments)
    {
        if(!mWriteFragments)
            return;

        fragments.forEach(x -> writeFragmentData(x, "UNFUSED"));
        fragments.forEach(x -> writeReadData(x.readId(), x.reads(), "UNFUSED"));
    }

    private void initialiseFragmentWriter()
    {
        if(!mWriteFragments)
            return;

        try
        {
            final String outputFileName = mConfig.formOutputFile("fusion_frags.csv");

            mFragmentWriter = createBufferedWriter(outputFileName, false);
            mFragmentWriter.write("ReadId,ReadCount,FusionGroup,Type,SameGeneSet,ScCount,HasSupp");

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final String prefix = se == SE_START ? "Start" : "End";
                mFragmentWriter.write(",Chr" + prefix);
                mFragmentWriter.write(",Orient" + prefix);
                mFragmentWriter.write(",JuncPos" + prefix);
                mFragmentWriter.write(",JuncOrient" + prefix);
                mFragmentWriter.write(",GeneSet" + prefix);
                mFragmentWriter.write(",Region" + prefix);
            }

            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric fragment data: {}", e.toString());
            return;
        }
    }

    public void writeIncompleteGroupReads(final List<FusionReadGroup> incompleteGroups)
    {
        incompleteGroups.forEach(x -> writeFusionReadData(x.ReadId, x.Reads, "INCOMPLETE_GROUPS"));
    }

    public synchronized void writeFragmentData(final FusionFragment fragment, final String fusionId)
    {
        if(!mWriteFragments)
            return;

        try
        {
            mFragmentWriter.write(String.format("%s,%d,%s,%s,%s,%d,%s",
                    fragment.readId(), fragment.reads().size(), fusionId, fragment.type(),
                    fragment.isSingleGeneCollection(),
                    fragment.reads().stream().filter(x -> x.SoftClipLengths[SE_START] > 0 || x.SoftClipLengths[SE_END] > 0).count(),
                    fragment.hasSuppAlignment()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                mFragmentWriter.write(String.format(",%s,%d,%d,%d,%d,%s",
                        fragment.chromosomes()[se], fragment.orientations()[se],
                        fragment.junctionPositions()[se], fragment.junctionOrientations()[se],
                        fragment.geneCollections()[se], fragment.regionMatchTypes()[se]));
            }

            mFragmentWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric fragment data: {}", e.toString());
        }
    }

    public synchronized void writeReadData(final String readId, final List<FusionRead> reads, final String groupStatus)
    {
        if(mWriteReads)
        {
            // not sure if will keep this
            writeFusionReadData(readId, reads, groupStatus);
        }
    }

    private void writeFusionReadData(final String readId, final List<FusionRead> reads, final String groupStatus)
    {
        if(mFusionReadWriter == null)
            return;

        try
        {
            for(final FusionRead read : reads)
            {
                mFusionReadWriter.write(String.format("%d,%s,%s,%s,%d,%d,%d,%s",
                        reads.size(), readId, groupStatus, read.Chromosome,
                        read.posStart(), read.posEnd(), read.Orientation, read.Cigar));

                /*
                mFusionReadWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%d,%s,%d",
                        read.isFirstOfPair(), read.HasSuppAlignment, read.isReadReversed(), read.isProperPair(),
                        read.SuppData != null, read.ReadBases, read.flags(), read.MateChromosome, read.MatePosStart));
                */

                mFusionReadWriter.write(String.format(",%d,%s,%s,%s,%s,%s,%d",
                        read.Flags, read.HasSuppAlignment, read.SuppData != null ? read.SuppData.asDelimStr() : "NONE",
                        read.BoundaryBases[SE_START], read.BoundaryBases[SE_END], read.MateChromosome, read.MatePosStart));

                mFusionReadWriter.write(String.format(",%d,%d,%s,%s,%s",
                        read.GeneCollections[SE_START], read.GeneCollections[SE_END],
                        read.IsGenicRegion[SE_START], read.IsGenicRegion[SE_END], read.HasInterGeneSplit));

                StringJoiner coordsStr = new StringJoiner(ITEM_DELIM);

                for(int[] coord : read.MappedCoords)
                {
                    coordsStr.add(String.format("%d:%d", coord[SE_START], coord[SE_END]));
                }

                mFusionReadWriter.write(String.format(",%s,%d,%d",
                        coordsStr, read.SoftClipLengths[SE_START], read.SoftClipLengths[SE_END]));

                // log the transcript exons affected, and the highest matching transcript
                StringJoiner transExonData = new StringJoiner(ITEM_DELIM);
                RegionMatchType topTransMatchType = read.getRegionMatchType(SE_START);

                if(topTransMatchType != NONE)
                {
                    for(final FusionTransExon transExonRef : read.getTransExonRefs(SE_START))
                    {
                        transExonData.add(String.format("%d:%d", transExonRef.TransId, transExonRef.ExonRank));
                    }
                }

                StringJoiner upperTransExonData = new StringJoiner(ITEM_DELIM);
                RegionMatchType upperTopTransMatchType = NONE;

                if(read.spansGeneCollections() && read.getTransExonRefs(SE_END) != null)
                {
                    upperTopTransMatchType = read.getRegionMatchType(SE_END);

                    for(FusionTransExon transExonRef : read.getTransExonRefs(SE_END))
                    {
                        upperTransExonData.add(String.format("%d:%d", transExonRef.TransId, transExonRef.ExonRank));
                    }
                }

                mFusionReadWriter.write(String.format(",%s,%s,%s,%s",
                        topTransMatchType, transExonData.toString().isEmpty() ? "NONE" : transExonData.toString(),
                        upperTopTransMatchType, upperTransExonData.toString().isEmpty() ? "NONE" : upperTransExonData.toString()));

                mFusionReadWriter.newLine();
            }

        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write chimeric read data: {}", e.toString());
            return;
        }
    }
}
