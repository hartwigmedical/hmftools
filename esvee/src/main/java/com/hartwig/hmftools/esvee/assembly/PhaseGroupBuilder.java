package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.SvConstants.PROXIMATE_DEL_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PROXIMATE_DUP_LENGTH;
import static com.hartwig.hmftools.esvee.common.AssemblySupport.hasMatchingFragment;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordant;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.DiscordantGroup;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.read.Read;

public class PhaseGroupBuilder
{
    private final SvConfig mConfig;
    private final Map<String, List<JunctionGroup>> mJunctionGroupMap;
    private final List<PhaseGroup> mPhaseGroups;

    public PhaseGroupBuilder(final SvConfig config, final Map<String, List<JunctionGroup>> junctionGroupMap)
    {
        mConfig = config;
        mJunctionGroupMap = junctionGroupMap;
        mPhaseGroups = Lists.newArrayList();

        // setting an index for each junction group allows easy access to an assembly's own group during look-ups
        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            for(int index = 0; index < junctionGroups.size(); ++index)
            {
                junctionGroups.get(index).setIndex(index);
            }
        }
    }

    private static final int LOG_COUNT = 10000;

    public void buildGroups()
    {
        int processed = 0;
        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            for(JunctionGroup junctionGroup : junctionGroups)
            {
                for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
                {
                    findLinkedAssemblies(junctionGroup, assembly);

                    ++processed;

                    if((processed % LOG_COUNT) == 0)
                    {
                        SV_LOGGER.debug("processed {} assemblies into {} phase groups",
                                processed, mPhaseGroups.size());
                    }
                }
            }
        }

        for(int i = 0; i < mPhaseGroups.size(); ++i)
        {
            mPhaseGroups.get(i).setId(i);
        }
    }

    public List<PhaseGroup> phaseGroups() { return mPhaseGroups; }

    private void findLinkedAssemblies(final JunctionGroup assemblyJunctionGroup, final JunctionAssembly assembly)
    {
        // for the given assembly, looks in all overlapping other junction groups (based on remote regions, ref-side soft-clips, and
        // the local junction group for indels) for other assemblies with shared reads
        if(assembly.remoteRegions().isEmpty() && assembly.refSideSoftClips().isEmpty() && !assembly.indel())
            return;

        PhaseGroup phaseGroup = assembly.phaseGroup(); // may have been set from an earlier assembly link
        boolean linksWithExisting = phaseGroup != null;

        Set<JunctionGroup> linkedJunctionGroups = Sets.newHashSet();

        for(RemoteRegion region : assembly.remoteRegions())
        {
            if(isFiltered(region))
                continue;

            List<JunctionGroup> overlappingJunctions = findOverlappingJunctionGroups(assemblyJunctionGroup, region);

            if(overlappingJunctions == null)
                continue;

            linkedJunctionGroups.addAll(overlappingJunctions);
        }

        if(!assembly.refSideSoftClips().isEmpty())
        {
            // if the assembly has candidate facing TI links, then add its own junction group - they will often share reads which are
            // discordant in one and a junction read in the other
            RefSideSoftClip refSideSoftClip = assembly.refSideSoftClips().get(0);

            if(positionWithin(refSideSoftClip.Position, assemblyJunctionGroup.minPosition(), assemblyJunctionGroup.maxPosition()))
            {
                linkedJunctionGroups.add(assemblyJunctionGroup);
            }
        }

        if(assembly.indel())
            linkedJunctionGroups.add(assemblyJunctionGroup);

        for(JunctionGroup junctionGroup : linkedJunctionGroups)
        {
            // the matching to other assemblies for each of this assembly's remote groups is purely for informational purposes
            // but in time may be replaced by actual linking to all expected remote mate reads

            for(JunctionAssembly otherAssembly : junctionGroup.junctionAssemblies())
            {
                if(assembly == otherAssembly)
                    continue;

                // for assemblies sharing read, the phase group scenarios are:
                // - assemblies are already in the same phase group
                // - no phase group exists for either so make a new one
                // - the other assembly already has an existing phase group and this one doesn't so just add it
                // - the other assembly already has an existing phase group, so does this one and so transfer this to the other
                if(phaseGroup != null && otherAssembly.phaseGroup() == phaseGroup)
                {
                    // already linked
                    continue;
                }

                boolean isLocalCandidateLink = (assemblyJunctionGroup == junctionGroup && isLocalAssemblyCandidate(assembly, otherAssembly));

                if(!isLocalCandidateLink && !assembliesShareReads(assembly, otherAssembly))
                    continue;

                // IDEA: first of all check the regions overlap, since for large junction groups there's a high chance they won't
                // and then all reads need to be checked in each
                // finds too many examples where they do share reads, maybe the 1K distance buffer is too small or they are sharing
                // supplementaries whose remote regions were purged

                /*
                boolean assemblyRegionsOverlap = assemblyJunctionGroup == junctionGroup || assembliesOverlap(assembly, otherAssembly);
                if(!assemblyRegionsOverlap)
                {
                    SV_LOGGER.debug("asm({}) and asm({}) share reads without any overlaps", assembly, otherAssembly);
                }
                */

                if(phaseGroup == null)
                {
                    if(otherAssembly.phaseGroup() != null)
                    {
                        phaseGroup = otherAssembly.phaseGroup();
                        phaseGroup.addAssembly(assembly);
                        linksWithExisting = true;
                    }
                    else
                    {
                        phaseGroup = new PhaseGroup(assembly, otherAssembly);
                    }
                }
                else
                {
                    if(otherAssembly.phaseGroup() != null)
                    {
                        // transfer to the other one, and do this for assemblies in this group
                        otherAssembly.phaseGroup().transferAssemblies(phaseGroup);
                        mPhaseGroups.remove(phaseGroup);
                        phaseGroup = otherAssembly.phaseGroup();
                        linksWithExisting = true;
                    }
                    else
                    {
                        phaseGroup.addAssembly(otherAssembly);
                    }
                }
            }
        }

        for(JunctionGroup junctionGroup : linkedJunctionGroups)
        {
            // the linked junction group may have discordant read groups which are required by this assembly
            for(DiscordantGroup discordantGroup : junctionGroup.discordantGroups())
            {
                // ignore any discordant group which actually shares the same reads
                if(discordantGroup.reads().stream().anyMatch(x -> assembly.hasReadSupport(x)))
                    continue;

                // as above, check that remote regions overlap before checking reads
                boolean matched = false;

                for(RemoteRegion remoteRegion : assembly.remoteRegions())
                {
                    if(remoteRegion.isSuppOnlyRegion())
                        continue;

                    if(!remoteRegion.overlaps(
                            discordantGroup.chromosome(), discordantGroup.minAlignedPosition(), discordantGroup.maxAlignedPosition(),
                            discordantGroup.orientation()))
                    {
                        continue;
                    }

                    if(remoteRegion.readIds().stream().anyMatch(x -> discordantGroup.hasFragment(x)))
                    {
                        matched = true;
                        break;
                    }
                }

                if(!matched)
                    continue;

                if(phaseGroup == null)
                    phaseGroup = new PhaseGroup(assembly, null);

                phaseGroup.addDiscordantGroup(discordantGroup);
            }
        }

        if(phaseGroup != null && !linksWithExisting)
        {
            mPhaseGroups.add(phaseGroup);
        }
    }

    private static boolean assembliesOverlap(final JunctionAssembly first, final JunctionAssembly second)
    {
        for(RemoteRegion region : first.remoteRegions())
        {
            int regionStart = region.start() - BAM_READ_JUNCTION_BUFFER;
            int regionEnd = region.start() - BAM_READ_JUNCTION_BUFFER;

            if(positionWithin(second.junction().Position, regionStart, regionEnd))
                return true;
        }

        for(RemoteRegion region : second.remoteRegions())
        {
            int regionStart = region.start() - BAM_READ_JUNCTION_BUFFER;
            int regionEnd = region.start() - BAM_READ_JUNCTION_BUFFER;

            if(positionWithin(first.junction().Position, regionStart, regionEnd))
                return true;
        }

        return false;
    }

    public static boolean isLocalAssemblyCandidate(final JunctionAssembly first, final JunctionAssembly second)
    {
        if(!first.junction().chromosome().equals(second.junction().chromosome()))
            return false;

        // assemblies must have DEL or DUP orientations, be within threshold distances of each other
        if(first.isForwardJunction() == second.isForwardJunction())
            return false;

        boolean firstIsLower = first.junction().Position <= second.junction().Position;
        boolean isDelType = firstIsLower == first.isForwardJunction();
        int junctionDistance = abs(first.junction().Position - second.junction().Position);

        if((isDelType && junctionDistance > PROXIMATE_DEL_LENGTH) || (!isDelType && junctionDistance > PROXIMATE_DUP_LENGTH))
            return false;

        // must have concordant reads with mates crossing the other junction
        JunctionAssembly lowerAssembly = firstIsLower ? first : second;
        JunctionAssembly upperAssembly = !firstIsLower ? first : second;
        Junction lowerJunction = firstIsLower ? first.junction() : second.junction();
        Junction upperJunction = !firstIsLower ? first.junction() : second.junction();

        if(lowerAssembly.support().stream().noneMatch(x -> isCrossingConcordantRead(x.read(), upperJunction, false)))
            return false;

        if(upperAssembly.support().stream().noneMatch(x -> isCrossingConcordantRead(x.read(), lowerJunction, true)))
            return false;

        return true;
    }

    private static boolean isCrossingConcordantRead(final Read read, final Junction junction, boolean requireLower)
    {
        if(isDiscordant(read) || read.isMateUnmapped() || !read.isPairedRead())
            return false;

        if(requireLower)
            return read.mateAlignmentEnd() < junction.Position;
        else
            return read.mateAlignmentStart() > junction.Position;
    }

    private static boolean assembliesShareReads(final JunctionAssembly first, final JunctionAssembly second)
    {
        // tests matching reads in both the junction reads and any extension reads (ie discordant)
        for(AssemblySupport support : first.support())
        {
            if(support.type() == JUNCTION_MATE)
                continue;

            if(hasMatchingFragment(second.support(), support.read()))
                return true;
        }

        for(AssemblySupport support : first.candidateSupport())
        {
            if(hasMatchingFragment(second.support(), support.read()))
                return true;
        }

        return false;
    }

    private List<JunctionGroup> findOverlappingJunctionGroups(final JunctionGroup assemblyJunctionGroup, final RemoteRegion region)
    {
        List<JunctionGroup> junctionGroups = mJunctionGroupMap.get(region.Chromosome);

        if(junctionGroups == null)
            return Collections.emptyList();

        int startIndex;

        if(assemblyJunctionGroup.chromosome().equals(region.Chromosome) && assemblyJunctionGroup.overlapsRemoteRegion(region))
        {
            startIndex = assemblyJunctionGroup.index();
        }
        else
        {
            // make use of the binary search for junction groups to find the closest index
            startIndex = JunctionGroup.binarySearch(region.start(), junctionGroups);
        }

        List<JunctionGroup> overlapGroups = Lists.newArrayList();

        for(int d = 0; d <= 1; ++d)
        {
            boolean searchUp = (d == 0);

            int index = searchUp ? startIndex + 1 : startIndex;

            while(index >= 0 && index < junctionGroups.size())
            {
                if(!junctionGroups.get(index).overlapsRemoteRegion(region))
                    break;

                overlapGroups.add(junctionGroups.get(index));

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }

        return overlapGroups;
    }

    private boolean isFiltered(final RemoteRegion region)
    {
        // ignore any remote region filtered out by config
        if(!mConfig.SpecificChrRegions.hasFilters())
            return false;

        if(mConfig.SpecificChrRegions.excludeChromosome(region.Chromosome))
            return true;

        if(mConfig.SpecificChrRegions.Regions.isEmpty())
            return false;

        return mConfig.SpecificChrRegions.Regions.stream().noneMatch(x -> x.overlaps(region));
    }
}
