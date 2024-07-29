package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsToStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.hasUnsetBases;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.CigarUtil;

public class AssemblyAlignment
{
    private final int mId;
    private final List<JunctionAssembly> mAssemblies;
    private final PhaseSet mPhaseSet;

    private final String mFullSequence;
    private final int mFullSequenceLength;

    private final Map<Integer,String> mSequenceOverlaps;
    private String mSequenceCigar;

    private final List<Breakend> mBreakends;

    private final Map<String,List<SupportRead>> mFragmentReadsMap;

    public AssemblyAlignment(final int id, final JunctionAssembly assembly) { this(id, assembly, null); }

    public AssemblyAlignment(final int id, final PhaseSet phaseSet) { this(id, null, phaseSet); }

    private AssemblyAlignment(final int id, final JunctionAssembly assembly, final PhaseSet phaseSet)
    {
        mId = id;

        if(assembly != null)
        {
            mAssemblies = List.of(assembly);
            mPhaseSet = null;
        }
        else
        {
            mPhaseSet = phaseSet;
            mAssemblies = mPhaseSet.assemblies();
        }

        mSequenceCigar = "";
        mSequenceOverlaps = Maps.newHashMap();
        mFragmentReadsMap = Maps.newHashMap();
        mBreakends = Lists.newArrayList();

        mAssemblies.forEach(x -> x.setAssemblyAlignmentInfo(info())); // set for output TSV only

        if(mAssemblies.stream().noneMatch(x -> hasUnsetBases(x)))
        {
            mFullSequence = buildSequenceData();
            mFragmentReadsMap.clear();
        }
        else
        {
            mFullSequence = null;
        }

        mFullSequenceLength = mFullSequence != null ? mFullSequence.length() : 0;
    }

    public int id() { return mId; }

    public boolean isValid() { return mFullSequence != null; }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }

    public PhaseSet phaseSet() { return mPhaseSet; }

    public List<Breakend> breakends() { return mBreakends; }
    public void addBreakend(final Breakend breakend) { mBreakends.add(breakend); }

    public String fullSequence() { return mFullSequence; }
    public int fullSequenceLength() { return mFullSequenceLength; }

    public String overlapBases(int sequenceIndex)
    {
        return mSequenceOverlaps.containsKey(sequenceIndex) ? mSequenceOverlaps.get(sequenceIndex) : "";
    }

    public String assemblyIds()
    {
        if(mAssemblies.size() == 1)
            return String.valueOf(mAssemblies.get(0).id());

        return mAssemblies.stream().map(x -> String.valueOf(x.id())).collect(Collectors.joining(ITEM_DELIM));
    }

    public String info()
    {
        if(mAssemblies.size() == 1)
            return mAssemblies.get(0).junction().coords();

        return mAssemblies.stream().map(x -> x.junction().coords()).collect(Collectors.joining("_"));
    }

    private String buildSingleAssemblySequenceData()
    {
        JunctionAssembly assembly = mAssemblies.get(0);

        int extensionLength = assembly.extensionLength();
        int refBaseLength = assembly.refBaseLength();

        if(assembly.isForwardJunction())
        {
            setAssemblyReadIndices(assembly, false, 0, refBaseLength);
        }
        else
        {
            setAssemblyReadIndices(assembly, false, extensionLength, refBaseLength);
        }

        if(assembly.isForwardJunction())
            mSequenceCigar = format("%dM%dS", refBaseLength, extensionLength);
        else
            mSequenceCigar = format("%dS%dM", extensionLength, refBaseLength);

        return assembly.formFullSequence();
    }

    private String buildSequenceData()
    {
        if(mAssemblies.size() == 1)
            return buildSingleAssemblySequenceData();

        // by convention if both ends of the phase set have the same orientation, then start with whichever is forward
        // and if they don't match then also start with the forward orientation
        // the exception being if both ends have outer SGLs (ie -ve to +ve orientation assemblies)

        List<JunctionAssembly> assemblies = mAssemblies;
        List<AssemblyLink> assemblyLinks = mPhaseSet.assemblyLinks();

        int assemblyCount = assemblies.size();
        int linkCount = assemblyLinks.size();

        // for a chain which starts with a facing link, reverse the order in which the sequence is added so as to start with ref bases
        if(assemblyLinks.get(0).type() == LinkType.FACING && assemblyLinks.get(linkCount - 1).type() != LinkType.FACING)
        {
            assemblies = Lists.newArrayList(assemblies);
            assemblyLinks = Lists.newArrayList(assemblyLinks);
            Collections.reverse(assemblies);
            Collections.reverse(assemblyLinks);
        }

        JunctionAssembly first = assemblies.get(0);
        JunctionAssembly last = assemblies.get(assemblyCount - 1);

        boolean startReversed, hasOuterExtensions;

        if(!first.isForwardJunction() && last.isForwardJunction() && mPhaseSet.hasFacingLinks())
        {
            startReversed = false;
            hasOuterExtensions = true;
        }
        else
        {
            startReversed = first.junction().isReverse();
            hasOuterExtensions = false;
        }

        SV_LOGGER.trace("building full alignment from {} assemblies, startRev({}) hasOuter({})",
                mAssemblies.size(), startReversed, hasOuterExtensions);

        StringBuilder fullSequence = new StringBuilder();

        boolean lastAddedReversed = false;
        int currentSeqLength = 0;
        boolean nextIsFacing = false;

        List<CigarElement> sequenceCigar = Lists.newArrayList();

        for(int i = 0; i < assemblyCount - 1; ++i)
        {
            AssemblyLink link = assemblyLinks.get(i);
            JunctionAssembly assembly = assemblies.get(i);

            boolean assemblyReversed;

            if(i == 0)
            {
                assemblyReversed = startReversed;

                if(hasOuterExtensions)
                {
                    // add on the extension sequence instead of the ref base sequence
                    String assemblyExtensionBases = assembly.formJunctionSequence();

                    fullSequence.append(assemblyExtensionBases);

                    // add the extra base since the junction index itself is not included in these extension bases
                    setAssemblyReadIndices(assembly, assemblyReversed, assemblyExtensionBases.length() + 1, 0);

                    logBuildInfo(assembly, currentSeqLength, assemblyExtensionBases.length(), assemblyReversed, "outer-ext-bases");

                    currentSeqLength = assemblyExtensionBases.length();

                    buildSequenceCigar(sequenceCigar, S, assemblyExtensionBases.length());
                }
                else
                {
                    String assemblyRefBases = startReversed ?
                            Nucleotides.reverseComplementBases(assembly.formRefBaseSequence()) : assembly.formRefBaseSequence();

                    fullSequence.append(assemblyRefBases);

                    setAssemblyReadIndices(assembly, assemblyReversed, 0, assemblyRefBases.length());

                    currentSeqLength += assemblyRefBases.length();

                    buildSequenceCigar(sequenceCigar, M, assemblyRefBases.length());

                    logBuildInfo(assembly, currentSeqLength, assemblyRefBases.length(), assemblyReversed, "ref-bases");
                }
            }
            else
            {
                if(nextIsFacing)
                {
                    nextIsFacing = false;

                    // ref bases for this segment have already been added so only set assembly indices
                    JunctionAssembly nextAssembly = assemblies.get(i + 1);

                    setAssemblyReadIndices(nextAssembly, lastAddedReversed, currentSeqLength, 0);
                    continue;
                }

                assemblyReversed = lastAddedReversed;
            }

            int overlapLength = link.overlapBases().length();
            String insertedBases = link.insertedBases();

            if(!insertedBases.isEmpty())
            {
                // keep inserted bases in the same direction as a full sequence
                fullSequence.append(assemblyReversed ? Nucleotides.reverseComplementBases(insertedBases) : insertedBases);

                currentSeqLength += insertedBases.length();

                buildSequenceCigar(sequenceCigar, I, insertedBases.length());
            }

            JunctionAssembly nextAssembly = assemblies.get(i + 1);
            boolean nextReversed = (assembly.junction().Orient == nextAssembly.junction().Orient && !assemblyReversed);

            String nextAssemblyRefBases = nextReversed ?
                    Nucleotides.reverseComplementBases(nextAssembly.formRefBaseSequence()) : nextAssembly.formRefBaseSequence();

            if(overlapLength >= nextAssemblyRefBases.length())
                return null;

            fullSequence.append(overlapLength > 0 ? nextAssemblyRefBases.substring(overlapLength) : nextAssemblyRefBases);

            currentSeqLength -= overlapLength;

            setAssemblyReadIndices(nextAssembly, nextReversed, currentSeqLength, nextAssemblyRefBases.length());

            logBuildInfo(nextAssembly, currentSeqLength, nextAssemblyRefBases.length(), assemblyReversed, "ref-bases");

            currentSeqLength += nextAssemblyRefBases.length();

            if(overlapLength > 0)
            {
                buildSequenceCigar(sequenceCigar, I, overlapLength);
                mSequenceOverlaps.put(currentSeqLength - 1, link.overlapBases());
            }

            buildSequenceCigar(sequenceCigar, M, nextAssembly.refBaseLength());

            if(hasOuterExtensions && i == assemblyCount - 2)
            {
                // add on the extension sequence for the last assembly
                String assemblyExtensionBases = nextAssembly.formJunctionSequence();

                fullSequence.append(assemblyExtensionBases);

                // setAssemblyReadIndices(assembly, assemblyReversed, currentSeqLength, 0);

                logBuildInfo(nextAssembly, currentSeqLength, assemblyExtensionBases.length(), assemblyReversed, "outer-ext-bases");

                currentSeqLength = assemblyExtensionBases.length();

                buildSequenceCigar(sequenceCigar, S, assemblyExtensionBases.length());
            }

            lastAddedReversed = nextReversed;

            nextIsFacing = true;
        }

        mSequenceCigar = cigarElementsToStr(sequenceCigar);

        return fullSequence.toString();
    }

    private static void logBuildInfo(
            final JunctionAssembly assembly, int currentSeqLength, int assemblyBaseLength, boolean isReversed, final String otherInfo)
    {
        if(SV_LOGGER.isTraceEnabled())
        {
            SV_LOGGER.trace("seqLength({} -> {}) adding assembly({}) {} {}",
                    currentSeqLength, currentSeqLength + assemblyBaseLength,
                    assembly.junction().coords(), isReversed ? "rev" : "fwd", otherInfo);
        }
    }

    private static void buildSequenceCigar(final List<CigarElement> elements, final CigarOperator operator, int length)
    {
        if(length == 0)
            return;

        if(elements.isEmpty())
        {
            elements.add(new CigarElement(length, operator));
            return;
        }

        int lastIndex = elements.size() - 1;
        CigarElement lastElement = elements.get(lastIndex);

        if(lastElement.getOperator() != operator)
        {
            elements.add(new CigarElement(length, operator));
            return;
        }

        elements.set(lastIndex, new CigarElement(lastElement.getLength() + length, operator));
    }

    private void setAssemblyReadIndices(
            final JunctionAssembly assembly, boolean isReversed, int existingSequenceLength, int assemblyRefBaseLength)
    {
        // for reads with forward orientation, their full sequence index is just the full sequence start offset + the read's junction distance
        // for reads in reversed assemblies, define their assembly index from the read's end index then going forward to the read start index

        // always start from the existing sequence index (ie existing length)
        // define existing length (E), assembly ref-base length (AR), read start distance from junction (JRSD), read assembly index (RI)
        // forward assemblies, not reversed: RI = E + AR - JRSD
        // reverse assemblies, not reversed: RI = E - JRSD (since is in relation to the junction being at the existing sequence end)
        // for reversed assemblies, the JRSD is negated ('JRSD) and the read's end is then used, so JRSD = -(JRSD) - read-length
        // forward assemblies, reversed: eg a +ve inversion has the second assembly reversed and added to the end
        //      RI = E - 'JRSD
        // reverse assemblies, reversed: eg a -ve inversion has the first assembly reversed and added to the start
        //      RI = E + AR - 'JRSD

        boolean includeAssemblyRefBaseLength = assembly.isForwardJunction() == !isReversed;

        int currentSeqEndIndex = existingSequenceLength - 1;

        if(includeAssemblyRefBaseLength)
            currentSeqEndIndex += assemblyRefBaseLength;

        for(SupportRead read : assembly.support())
        {
            int juncReadStartDistance = read.junctionReadStartDistance();

            Orientation fullSeqOrientation = isReversed ? Orientation.REVERSE : Orientation.FORWARD;

            int fullSeqIndex = currentSeqEndIndex;

            if(!isReversed)
            {
                // typical scenario - eg for a junction read, this will bring the read's start index to earlier (lower) than the assembly index
                fullSeqIndex -= juncReadStartDistance;
            }
            else
            {
                // eg a junction read with JRSD = 100, ie the read starting 100 bases before the junction, read end at +50 into extension bases
                // after reversing, the read starts (from its end at -50 relative to this assembly's reversed junction end, extending to +100 past it
                fullSeqIndex += juncReadStartDistance - read.baseLength() + 1;
            }

            SupportRead matchedRead = findMatchingSupportRead(read);

            if(matchedRead != null)
            {
                // always use any earlier orientation
                fullSeqOrientation = matchedRead.fullAssemblyOrientation();

                // use junction reads over discordant reads, use earlier junction reads over later, and use the discordant
                // read closest to a junction
                boolean useMatchedReadInfo;

                if(read.type() == matchedRead.type() && !read.type().isSplitSupport())
                {
                    useMatchedReadInfo = juncReadStartDistance > matchedRead.junctionReadStartDistance();
                }
                else
                {
                    useMatchedReadInfo = true;
                }

                if(useMatchedReadInfo)
                {
                    fullSeqIndex = matchedRead.fullAssemblyIndexStart();
                }
                else
                {
                    matchedRead.setFullAssemblyInfo(fullSeqIndex, fullSeqOrientation);
                }
            }

            read.setFullAssemblyInfo(fullSeqIndex, fullSeqOrientation);
        }

        // TODO: any secondaries need their support reads' full sequence assembly indices to be set too

    }

    private SupportRead findMatchingSupportRead(final SupportRead read)
    {
        List<SupportRead> fragmentReads = mFragmentReadsMap.get(read.id());

        if(fragmentReads == null)
        {
            mFragmentReadsMap.put(read.id(), Lists.newArrayList(read));
            return null;
        }

        SupportRead matchedRead = fragmentReads.stream().filter(x -> x.firstInPair() == read.firstInPair()).findFirst().orElse(null);
        fragmentReads.add(read);
        return matchedRead;
    }

    public String assemblyCigar() { return mSequenceCigar; }

    public String toString()
    {
        return format("%s length(%d) breakends(%d)", info(), mFullSequenceLength, mBreakends.size());
    }
}
