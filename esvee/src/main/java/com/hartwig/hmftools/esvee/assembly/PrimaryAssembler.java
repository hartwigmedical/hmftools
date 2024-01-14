package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_MISMATCH_READS;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_MISMATCH_TOTAL_QUAL;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.buildFromJunctionReads;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.purgeLowSupport;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.common.AssemblyMismatchSplitter;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.output.ResultsWriter;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadFilters;

public class PrimaryAssembler
{
    private final SvConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final Junction mJunction;

    private int mNextAssemblyNumber = 1;

    public PrimaryAssembler(final SvConfig config, final ResultsWriter resultsWriter, final Junction junction)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mJunction = junction;
    }

    public List<JunctionAssembly> processJunction(final List<Read> rawReads)
    {
        // FIXME:
        final List<Read> realignedReads = rawReads.stream()
                // .map(alignment -> realignForJunction(alignment, mJunction))
                .collect(Collectors.toList());

        final List<Read> withLowQAlignments = realignedReads.stream()
                .filter(alignment -> ReadFilters.recordSoftClipsNearJunction(alignment, mJunction)) // mCounters.ReadsSoftClippedAtJunction
                .collect(Collectors.toList());

        final List<Read> filteredAlignments = withLowQAlignments.stream()
                .filter(alignment -> ReadFilters.isRecordAverageQualityPastJunctionAbove(alignment, mJunction, SvConstants.AVG_BASE_QUAL_THRESHOLD)) // mCounters.ReadsPassingJunctionQualityThreshold
                .filter(alignment -> ReadFilters.hasAcceptableMapQ(alignment, SvConstants.READ_FILTER_MIN_JUNCTION_MAPQ)) // mCounters.HasAcceptableMapQ
                .filter(ReadFilters::isNotBadlyMapped) // mCounters.WellMapped
                .collect(Collectors.toList());

        if(filteredAlignments.size() < PRIMARY_ASSEMBLY_MIN_MISMATCH_READS)
            return List.of();

        List<JunctionAssembly> initialAssemblies = createInitialAssemblies(filteredAlignments);

        initialAssemblies.forEach(x -> mResultsWriter.writeAssembly(x));

        return initialAssemblies;
    }

    private String nextAssemblyName()
    {
        // consider naming based on initial length and support? try to make typically deterministic
        return String.format("%s:%s%s:%s", mJunction.Chromosome, mJunction.Position,
                mJunction.direction() == Direction.FORWARDS ? "F" : "R", mNextAssemblyNumber++);
    }

    private List<JunctionAssembly> createInitialAssemblies(final List<Read> junctionReads)
    {
        JunctionAssembly junctionSequence = buildFromJunctionReads(mJunction, junctionReads, true);

        if(junctionSequence.length() < PRIMARY_ASSEMBLY_MIN_LENGTH)
            return Collections.emptyList();

        // filter
        boolean hasValidMismatches = purgeLowSupport(
                junctionSequence, PRIMARY_ASSEMBLY_MIN_MISMATCH_READS, PRIMARY_ASSEMBLY_MIN_MISMATCH_TOTAL_QUAL);

        List<JunctionAssembly> junctionSequences;

        if(hasValidMismatches)
        {
            AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
            junctionSequences = splitter.splitOnMismatches(PRIMARY_ASSEMBLY_MIN_LENGTH);
        }
        else
        {
            junctionSequences = List.of(junctionSequence);
        }

        // extend these sequences in the direction away from the junction
        junctionSequences.forEach(x -> x.expandReferenceBases());



        return junctionSequences;
    }

}
