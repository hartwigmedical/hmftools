package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_MISMATCH_READS;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_READ_MAX_BASE_MISMATCH;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.buildFromJunctionReads;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
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
        List<Read> junctionReads = Lists.newArrayList();
        List<Read> otherReads = Lists.newArrayList();

        for(Read read : rawReads)
        {
            boolean isJunctionRead = true;

            try
            {
                if(!ReadFilters.recordSoftClipsNearJunction(read, mJunction))
                {
                    isJunctionRead = false;
                }
                else if(!ReadFilters.isRecordAverageQualityPastJunctionAbove(read, mJunction, SvConstants.AVG_BASE_QUAL_THRESHOLD))
                {
                    isJunctionRead = false;
                }
                else if(!ReadFilters.hasAcceptableMapQ(read, SvConstants.READ_FILTER_MIN_JUNCTION_MAPQ))
                {
                    isJunctionRead = false;
                }
                else if(ReadFilters.isBadlyMapped(read))
                {
                    isJunctionRead = false;
                }

                if(isJunctionRead)
                    junctionReads.add(read);
                else
                    otherReads.add(read);
            }
            catch(Exception e)
            {
                SV_LOGGER.error("error filtering read: {}", read.toString());
            }
        }

        /*
            // FIXME: this converts indels near this specific junction to soft-clips. Does it really make a difference
            // given that this was already checked for all reads? When would it still be required and be beneficial?
            // NOTE: has no mention in doco

        final List<Read> realignedReads = rawReads.stream()
                //.map(alignment -> realignForJunction(alignment, mJunction))
                .collect(Collectors.toList());
        final List<Read> withLowQAlignments = rawReads.stream()
                .filter(alignment -> ReadFilters.recordSoftClipsNearJunction(alignment, mJunction)) // mCounters.ReadsSoftClippedAtJunction
                .collect(Collectors.toList());

        final List<Read> filteredAlignments = withLowQAlignments.stream()
                .filter(alignment -> ReadFilters.isRecordAverageQualityPastJunctionAbove(alignment, mJunction, SvConstants.AVG_BASE_QUAL_THRESHOLD)) // mCounters.ReadsPassingJunctionQualityThreshold
                .filter(alignment -> ReadFilters.hasAcceptableMapQ(alignment, SvConstants.READ_FILTER_MIN_JUNCTION_MAPQ)) // mCounters.HasAcceptableMapQ
                .filter(ReadFilters::isNotBadlyMapped) // mCounters.WellMapped
                .collect(Collectors.toList());
        */

        if(junctionReads.size() < PRIMARY_ASSEMBLY_MIN_MISMATCH_READS)
            return List.of();

        List<JunctionAssembly> initialAssemblies = createInitialAssemblies(junctionReads);

        // look for further support for non-junction reads
        for(Read read : otherReads)
        {
            for(JunctionAssembly assembly : initialAssemblies)
            {
                if(assembly.checkReadMatches(read, PRIMARY_ASSEMBLY_READ_MAX_BASE_MISMATCH))
                {
                    assembly.addRead(read, false);
                }
            }
        }

        // dedup assemblies for this junction based on overlapping read support




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

        // no filtering of the initial sequence and instead rely on the sequence splitting to do this with all initial mismatches preserved
        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
        List<JunctionAssembly> junctionSequences = splitter.splitOnMismatches(PRIMARY_ASSEMBLY_MIN_LENGTH);

        // extend these sequences in the direction away from the junction
        junctionSequences.forEach(x -> x.expandReferenceBases());

        return junctionSequences;
    }

}
