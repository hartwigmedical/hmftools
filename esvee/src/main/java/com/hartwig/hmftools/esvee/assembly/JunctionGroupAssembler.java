package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.processor.PrimaryAssemblyResult;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;

import htsjdk.samtools.SAMRecord;

public class JunctionGroupAssembler
{
    private final SvConfig mConfig;
    private final JunctionGroup mJunctionGroup;
    private final BamReader mBamReader;

    private final Map<String,Read> mReadGroupMap;

    private final List<PrimaryAssemblyResult> mPrimaryAssemblyResults;

    public JunctionGroupAssembler(final SvConfig config, final BamReader bamReader, final JunctionGroup junctionGroup)
    {
        mConfig = config;
        mJunctionGroup = junctionGroup;
        mBamReader = bamReader;

        mReadGroupMap = Maps.newHashMap();
        mPrimaryAssemblyResults = Lists.newArrayList();
    }

    public List<PrimaryAssemblyResult> primaryAssemblyResults() { return mPrimaryAssemblyResults; }

    public void clear()
    {
        mJunctionGroup.clearCandidateReads();
        mPrimaryAssemblyResults.clear();
        mReadGroupMap.clear();
    }

    public void run()
    {
        int sliceStart = mJunctionGroup.minPosition() - BAM_READ_JUNCTION_BUFFER;
        int sliceEnd = mJunctionGroup.maxPosition() + BAM_READ_JUNCTION_BUFFER;

        SV_LOGGER.trace("junctionGroup({}:{}-{} count={}) slice",
                mJunctionGroup.chromosome(), sliceStart, sliceEnd, mJunctionGroup.count());

        mBamReader.sliceBam(mJunctionGroup.chromosome(), sliceStart, sliceEnd, this::processRecord);

        SV_LOGGER.debug("junctionGroup({}:{}-{} count={}) slice complete, readCount({})",
                mJunctionGroup.chromosome(), sliceStart, sliceEnd, mJunctionGroup.count(), mJunctionGroup.candidateReadCount());

        // now pass applicable reads to each primary assembler
        for(Junction junction : mJunctionGroup.junctions())
        {
            PrimaryAssembler primaryAssembler = new PrimaryAssembler(mConfig, junction);

            List<Read> junctionCandidateReads = mJunctionGroup.candidateReads().stream()
                    .filter(x -> AlignmentFilters.alignmentCrossesJunction(x, junction))
                    .collect(Collectors.toList());

            if(junctionCandidateReads.isEmpty())
                continue;

            List<PrimaryAssembly> candidateAssemblies = primaryAssembler.processJunction(junctionCandidateReads);

            mPrimaryAssemblyResults.add(new PrimaryAssemblyResult(junction, primaryAssembler.getCounters(), candidateAssemblies));
        }
    }

    private static final int MATE_READ_BUFFER = 200;

    private void processRecord(final SAMRecord record)
    {
        // CHECK: do in SvPrep if worthwhile
        if(!AlignmentFilters.isRecordAverageQualityAbove(record.getBaseQualities(), SvConstants.AVG_BASE_QUAL_THRESHOLD))
            return;

        // mReadRescue::rescueRead) // CHECK: not required

        // mNormaliser::normalise() on each record

        Read read = new Read(record);

        mJunctionGroup.addCandidateRead(read);

        // link first and second in pair if within the same group
        if(read.isMateMapped() && read.mateChromosome().equals(read.getChromosome())
        && positionWithin(read.mateAlignmentStart(), mJunctionGroup.minPosition() - MATE_READ_BUFFER, mJunctionGroup.maxPosition()))
        {
            Read mateRead = mReadGroupMap.get(read.getName());

            if(mateRead != null)
            {
                mReadGroupMap.remove(read.getName());
                mateRead.setMateRead(read);
                read.setMateRead(mateRead);
            }
            else
            {
                mReadGroupMap.put(read.getName(), read);
            }
        }
    }
}
