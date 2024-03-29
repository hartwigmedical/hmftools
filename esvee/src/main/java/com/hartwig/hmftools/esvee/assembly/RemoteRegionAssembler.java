package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.common.CommonUtils.deleteInterimFile;
import static com.hartwig.hmftools.esvee.common.CommonUtils.initialise;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordantFragment;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.types.AssemblyLink;
import com.hartwig.hmftools.esvee.types.AssemblySupport;
import com.hartwig.hmftools.esvee.types.Junction;
import com.hartwig.hmftools.esvee.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.types.RemoteRegion;
import com.hartwig.hmftools.esvee.types.SupportType;

import htsjdk.samtools.SAMRecord;

public class RemoteRegionAssembler
{
    private final AssemblyConfig mConfig;
    private final BamReader mBamReader;

    private ChrBaseRegion mRemoteRegion;
    private final Map<String,Read> mSourceReads;
    private final List<Read> mMatchedRemoteReads;

    public RemoteRegionAssembler(final AssemblyConfig config, final BamReader bamReader)
    {
        mConfig = config;
        mBamReader = bamReader;

        mRemoteRegion = null;
        mSourceReads = Maps.newHashMap();
        mMatchedRemoteReads = Lists.newArrayList();
    }

    public static boolean isExtensionCandidateAssembly(final JunctionAssembly assembly)
    {
        if(assembly.refBaseTrimLength() < MIN_VARIANT_LENGTH)
            return false;

        int maxExtBaseMatchCount = 0;
        int secondExtBaseMatchCount = 0;
        int remoteJuncMates = 0;

        for(AssemblySupport support : assembly.support())
        {
            if(support.read().isReference()) // for now no reference support
                return false;

            Read read = support.read();

            if(support.type() != SupportType.JUNCTION)
                continue;

            if(support.junctionMatches() > maxExtBaseMatchCount)
            {
                secondExtBaseMatchCount = maxExtBaseMatchCount; // promote the second highest
                maxExtBaseMatchCount = support.junctionMatches();
            }
            else if(support.junctionMatches() > secondExtBaseMatchCount)
            {
                secondExtBaseMatchCount = support.junctionMatches();
            }

            if(read.isSupplementary() || read.isMateUnmapped())
                continue;

            boolean matePastJunction = (read.orientation() == POS_ORIENT) == assembly.isForwardJunction();

            if(isDiscordantFragment(read))
            {
                if(matePastJunction)
                    ++remoteJuncMates;
            }
            else
            {
                if((assembly.isForwardJunction() && read.mateAlignmentStart() > assembly.junction().Position)
                || (!assembly.isForwardJunction() && read.mateAlignmentEnd() < assembly.junction().Position))
                {
                    return false;
                }
            }
        }

        if(secondExtBaseMatchCount < MIN_VARIANT_LENGTH)
            return false;

        if(remoteJuncMates < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return false;

        return true;
    }

    public static boolean assemblyOverlapsRemoteRegion(final JunctionAssembly assembly, final RemoteRegion remoteRegion)
    {
        return positionsOverlap(assembly.minAlignedPosition(), assembly.maxAlignedPosition(), remoteRegion.start(), remoteRegion.end());
    }

    public void formRemoteAssembly(final JunctionAssembly assembly, final RemoteRegion remoteRegion, final List<Read> sourceReads)
    {
        mRemoteRegion = remoteRegion;

        SV_LOGGER.trace("remote region({}) slice", mRemoteRegion);

        mMatchedRemoteReads.clear();
        sourceReads.forEach(x -> mSourceReads.put(x.id(), x));

        mBamReader.sliceBam(mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end(), this::processRecord);

        SV_LOGGER.debug("remote region({}) sourcedReads(matched={} unmatched={})",
                mRemoteRegion, mMatchedRemoteReads.size(), mSourceReads.size());

        if(mMatchedRemoteReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return;

        // form a remote ref-base assembly from these reads
        int remoteRegionStart = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentStart()).min().orElse(0);
        int remoteRegionEnd = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentEnd()).max().orElse(0);

        byte[] refGenomeBases = mConfig.RefGenome.getBases(remoteRegion.Chromosome, remoteRegionStart, remoteRegionEnd);

        int remoteApproxJunctionPosition = remoteRegion.isForward() ? remoteRegionEnd : remoteRegionStart;
        Junction remoteJunction = new Junction(remoteRegion.Chromosome, remoteApproxJunctionPosition, remoteRegion.orientation());

        byte[] refBaseQuals = new byte[refGenomeBases.length];
        initialise(refBaseQuals, (byte)(LOW_BASE_QUAL_THRESHOLD + 1));

        JunctionAssembly remoteAssembly = new JunctionAssembly(
                remoteJunction, refGenomeBases, refBaseQuals, Collections.emptyList(), Collections.emptyList());

        AssemblyLinker assemblyLinker = new AssemblyLinker();

        AssemblyLink assemblyLink = assemblyLinker.tryAssemblyOverlap(assembly, remoteAssembly);

        if(assemblyLink != null)
        {
            SV_LOGGER.debug("assembly({}) links with remote region({}) matchedReads({})",
                    assembly, remoteRegion.toString(), mMatchedRemoteReads.size());
        }
    }

    private void processRecord(final SAMRecord record)
    {
        Read sourceRead = mSourceReads.remove(record.getReadName());

        if(sourceRead == null)
            return;

        Read remoteRead = new Read(record);

        if(mBamReader.currentIsReferenceSample())
            remoteRead.markReference();

        mMatchedRemoteReads.add(remoteRead);
    }
}
