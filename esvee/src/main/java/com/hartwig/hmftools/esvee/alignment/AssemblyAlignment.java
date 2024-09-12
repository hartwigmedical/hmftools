package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsToStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.extractInsertSequence;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.hasUnsetBases;
import static com.hartwig.hmftools.esvee.assembly.types.LinkType.FACING;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.AssemblyUtils;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class AssemblyAlignment
{
    private int mId;
    private final List<JunctionAssembly> mAssemblies;
    private final PhaseSet mPhaseSet;

    private String mFullSequence;
    private int mFullSequenceLength;

    private final Map<Integer,String> mSequenceOverlaps;
    private String mSequenceCigar;

    private final List<Breakend> mBreakends;
    private final List<Integer> mLinkIndices; // indices within the full sequence of each break junction

    private final Map<String,List<SupportRead>> mFragmentReadsMap;

    public AssemblyAlignment(final JunctionAssembly assembly) { this(assembly, null); }

    public AssemblyAlignment(final PhaseSet phaseSet) { this(null, phaseSet); }

    private AssemblyAlignment(final JunctionAssembly assembly, final PhaseSet phaseSet)
    {
        mId = -1;

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
        mLinkIndices = Lists.newArrayList();

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

    public void setId(int id) { mId = id; }
    public int id() { return mId; }

    public boolean isValid() { return mFullSequence != null && mFullSequenceLength < 100000; }

    public List<JunctionAssembly> assemblies() { return mAssemblies; }

    public PhaseSet phaseSet() { return mPhaseSet; }
    public boolean isMerged() { return mPhaseSet != null && mPhaseSet.merged(); }

    public List<Breakend> breakends() { return mBreakends; }
    public void addBreakend(final Breakend breakend) { mBreakends.add(breakend); }

    public String fullSequence() { return mFullSequence; }
    public int fullSequenceLength() { return mFullSequenceLength; }

    public String assemblyCigar() { return mSequenceCigar; }
    public List<Integer> linkIndices() { return mLinkIndices; }

    public Map<Integer,String> sequenceOverlaps() { return mSequenceOverlaps; }

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
            setAssemblyReadIndices(assembly, false, refBaseLength - 1);
        }
        else
        {
            setAssemblyReadIndices(assembly, false, extensionLength);
        }

        if(assembly.isForwardJunction())
        {
            mSequenceCigar = format("%dM%dS", refBaseLength, extensionLength);
            mLinkIndices.add(refBaseLength);
        }
        else
        {
            mSequenceCigar = format("%dS%dM", extensionLength, refBaseLength);
            mLinkIndices.add(extensionLength);
        }

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
        boolean hasFacingAtStart = false;
        boolean hasFacingAtEnd = false;

        if(mPhaseSet.hasFacingLinks())
        {
            hasFacingAtStart = assemblyLinks.get(0).type() == FACING;
            hasFacingAtEnd = assemblyLinks.get(linkCount - 1).type() == FACING;
        }

        boolean hasOuterExtensions = hasFacingAtStart && hasFacingAtEnd;

        boolean reverseLinks = hasFacingAtStart && !hasFacingAtEnd;

        if(reverseLinks)
        {
            hasFacingAtStart = false;
            hasFacingAtEnd = true;
            assemblies = Lists.newArrayList(assemblies);
            assemblyLinks = Lists.newArrayList(assemblyLinks);
            Collections.reverse(assemblies);
            Collections.reverse(assemblyLinks);
        }

        boolean startReversed = hasOuterExtensions ? false : assemblies.get(0).junction().isReverse();

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

                if(hasFacingAtStart)
                {
                    // add on the extension sequence instead of the ref base sequence
                    String assemblyExtensionBases = assembly.formJunctionSequence();
                    int assemblyExtensionBaseLength = assemblyExtensionBases.length();

                    fullSequence.append(assemblyExtensionBases);

                    // add the extra base since the junction index itself is not included in these extension bases
                    setAssemblyReadIndices(assembly, assemblyReversed, assemblyExtensionBaseLength);

                    logBuildInfo(assembly, currentSeqLength, assemblyExtensionBaseLength, assemblyReversed, "outer-ext-bases");

                    currentSeqLength = assemblyExtensionBaseLength;

                    buildSequenceCigar(sequenceCigar, S, assemblyExtensionBaseLength);
                    mLinkIndices.add(assemblyExtensionBaseLength);
                }
                else
                {
                    String assemblyRefBases = startReversed ?
                            Nucleotides.reverseComplementBases(assembly.formRefBaseSequence()) : assembly.formRefBaseSequence();
                    int assemblyRefBaseLength = assemblyRefBases.length();

                    fullSequence.append(assemblyRefBases);

                    setAssemblyReadIndices(assembly, assemblyReversed, assemblyRefBaseLength - 1);

                    currentSeqLength += assemblyRefBaseLength;

                    buildSequenceCigar(sequenceCigar, M, assemblyRefBaseLength);

                    logBuildInfo(assembly, currentSeqLength, assemblyRefBaseLength, assemblyReversed, "ref-bases");
                }
            }
            else
            {
                if(nextIsFacing)
                {
                    nextIsFacing = false;

                    // ref bases for this segment have already been added so only set assembly indices
                    JunctionAssembly nextAssembly = assemblies.get(i + 1);

                    setAssemblyReadIndices(nextAssembly, lastAddedReversed, currentSeqLength);

                    if(i == assemblyCount - 2)
                    {
                        // add on the extension sequence for the last assembly
                        addFinalFacingExtensionBases(
                                nextAssembly, fullSequence, sequenceCigar, currentSeqLength, lastAddedReversed, assemblyLinks);

                        currentSeqLength = nextAssembly.extensionLength();
                    }

                    continue;
                }

                assemblyReversed = lastAddedReversed;
            }

            int overlapLength = link.overlapBases().length();
            int insertedBaseLength = link.insertedBases().length();

            JunctionAssembly nextAssembly = assemblies.get(i + 1);

            // typically reverse any assemblies which come in on the +ve side, unless this is a facing link (from outer extensions)
            boolean nextReversed = nextAssembly.isForwardJunction() && link.type() != FACING;

            if(insertedBaseLength > 0)
            {
                // keep inserted bases in the same direction as the assembly is added
                String insertedBases = extractInsertSequence(assembly, assemblyReversed, nextAssembly, nextReversed, insertedBaseLength);

                fullSequence.append(insertedBases);

                currentSeqLength += insertedBaseLength;

                buildSequenceCigar(sequenceCigar, I, insertedBaseLength);
            }

            String nextAssemblyRefBases = nextReversed ?
                    Nucleotides.reverseComplementBases(nextAssembly.formRefBaseSequence()) : nextAssembly.formRefBaseSequence();

            if(overlapLength >= nextAssemblyRefBases.length())
                return null;

            if(overlapLength > 0) // remove the repeated / overlapping bases from the next ref bases
                nextAssemblyRefBases = nextAssemblyRefBases.substring(overlapLength);

            fullSequence.append(nextAssemblyRefBases);
            mLinkIndices.add(currentSeqLength);

            int nextAssemblyJunctionIndex = link.type() != FACING ? currentSeqLength : currentSeqLength + nextAssemblyRefBases.length();
            setAssemblyReadIndices(nextAssembly, nextReversed, nextAssemblyJunctionIndex);

            logBuildInfo(nextAssembly, currentSeqLength, nextAssemblyRefBases.length(), assemblyReversed, "ref-bases");

            if(overlapLength > 0)
            {
                buildSequenceCigar(sequenceCigar, I, overlapLength);
                mSequenceOverlaps.put(currentSeqLength - 1, link.overlapBases());
            }

            currentSeqLength += nextAssemblyRefBases.length(); // this will have been reduced by any overlapping bases already

            buildSequenceCigar(sequenceCigar, M, nextAssembly.refBaseLength());

            if(hasFacingAtEnd && i == assemblyCount - 2)
            {
                // add on the extension sequence for the last assembly
                addFinalFacingExtensionBases(nextAssembly, fullSequence, sequenceCigar, currentSeqLength, assemblyReversed, assemblyLinks);
                currentSeqLength = nextAssembly.extensionLength();
            }

            lastAddedReversed = nextReversed;

            nextIsFacing = true;
        }

        mSequenceCigar = cigarElementsToStr(sequenceCigar);

        return fullSequence.toString();
    }

    private void addFinalFacingExtensionBases(
            final JunctionAssembly assembly, final StringBuilder fullSequence, final List<CigarElement> sequenceCigar,
            int currentSeqLength, boolean assemblyReversed, final List<AssemblyLink> assemblyLinks)
    {
        // if the assembly is linked to another assembly at the start then repeat those other ref bases, otherwise take the extension
        String assemblyExtensionBases;
        CigarOperator extensionCigarType;
        String extensionInfo;

        JunctionAssembly matchedAssembly = assemblyLinks.get(0).findMatchedAssembly(assembly);

        if(matchedAssembly != null)
        {
            /* omit the repeated ref bases to avoid non-identical breakends being called

            JunctionAssembly linkedAssembly = assemblyLinks.get(0).otherAssembly(matchedAssembly);
            boolean nextReversed = linkedAssembly.isForwardJunction();

            assemblyExtensionBases = nextReversed ?
                    Nucleotides.reverseComplementBases(linkedAssembly.formRefBaseSequence()) : linkedAssembly.formRefBaseSequence();

            extensionCigarType = M;
            extensionInfo = "linked ref-bases";
            */
            return;
        }
        else
        {
            assemblyExtensionBases = assembly.isReverseJunction() ?
                    Nucleotides.reverseComplementBases(assembly.formJunctionSequence()) : assembly.formJunctionSequence();

            extensionCigarType = S;
            extensionInfo = "outer-ext-bases";
        }

        fullSequence.append(assemblyExtensionBases);

        logBuildInfo(assembly, currentSeqLength, assemblyExtensionBases.length(), assemblyReversed, extensionInfo);

        buildSequenceCigar(sequenceCigar, extensionCigarType, assemblyExtensionBases.length());
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

    private void setAssemblyReadIndices(final JunctionAssembly assembly, boolean isReversed, int assemblyJunctionSeqIndex)
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

        int currentSeqEndIndex = assemblyJunctionSeqIndex;

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

    public void updateSequenceInfo(final String newSequence, final Map<Integer,String> newSequenceOverlaps, int primaryOffsetAdjust)
    {
        mFullSequence = newSequence;
        mFullSequenceLength = newSequence.length();

        for(Map.Entry<Integer,String> entry : mSequenceOverlaps.entrySet())
        {
            entry.setValue(entry.getValue() + primaryOffsetAdjust);
        }

        mSequenceOverlaps.putAll(newSequenceOverlaps);
    }

    public String toString()
    {
        return format("%s length(%d) breakends(%d)", info(), mFullSequenceLength, mBreakends.size());
    }
}
