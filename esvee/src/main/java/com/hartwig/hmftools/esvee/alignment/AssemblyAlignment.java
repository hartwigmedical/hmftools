package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

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

        mFullSequence = buildSequenceData();
        mFullSequenceLength = mFullSequence != null ? mFullSequence.length() : 0;

        mBreakends = Lists.newArrayList();

        mAssemblies.forEach(x -> x.setAssemblyAlignmentInfo(info())); // set for output TSV only
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

        int assemblyCount = mAssemblies.size();
        JunctionAssembly first = mAssemblies.get(0);
        JunctionAssembly last = mAssemblies.get(assemblyCount - 1);

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

        StringBuilder fullSequence = new StringBuilder();
        StringBuilder sequenceCigar = new StringBuilder();

        boolean lastAddedReversed = false;
        int currentFullSeqLength = 0;

        for(int i = 0; i < assemblyCount - 1; i += 2) // to skip past facing links
        {
            AssemblyLink link = mPhaseSet.assemblyLinks().get(i);
            JunctionAssembly assembly = mPhaseSet.assemblies().get(i);

            boolean assemblyReversed;

            if(i == 0)
            {
                assemblyReversed = startReversed;

                if(hasOuterExtensions)
                {
                    // add on the extension sequence instead of the ref base sequence
                    String assemblyExtensionBases = assembly.formJunctionSequence();

                    fullSequence.append(assemblyExtensionBases);

                    setAssemblyReadIndices(assembly, assemblyReversed, assemblyExtensionBases.length(), 0);

                    currentFullSeqLength = assemblyExtensionBases.length();

                    sequenceCigar.append(format("%dS", assemblyExtensionBases.length()));
                }
                else
                {
                    String assemblyRefBases = startReversed ?
                            Nucleotides.reverseComplementBases(assembly.formRefBaseSequence()) : assembly.formRefBaseSequence();

                    fullSequence.append(assemblyRefBases);

                    setAssemblyReadIndices(assembly, assemblyReversed, 0, assemblyRefBases.length());

                    currentFullSeqLength = assemblyRefBases.length();

                    sequenceCigar.append(format("%dM", assemblyRefBases.length()));
                }
            }
            else
            {
                assemblyReversed = lastAddedReversed;
            }

            int overlapLength = link.overlapBases().length();
            String insertedBases = link.insertedBases();

            if(!insertedBases.isEmpty())
            {
                // keep inserted bases in the same direction as a full sequence
                fullSequence.append(assemblyReversed ? Nucleotides.reverseComplementBases(insertedBases) : insertedBases);

                currentFullSeqLength += insertedBases.length();

                sequenceCigar.append(format("%dI", insertedBases.length()));
            }

            JunctionAssembly nextAssembly = mPhaseSet.assemblies().get(i + 1);
            boolean nextReversed = (assembly.junction().Orient == nextAssembly.junction().Orient && !assemblyReversed);

            String nextAssemblyRefBases = nextReversed ?
                    Nucleotides.reverseComplementBases(nextAssembly.formRefBaseSequence()) : nextAssembly.formRefBaseSequence();

            if(overlapLength >= nextAssemblyRefBases.length())
                return null;

            fullSequence.append(overlapLength > 0 ? nextAssemblyRefBases.substring(overlapLength) : nextAssemblyRefBases);

            currentFullSeqLength -= overlapLength;

            setAssemblyReadIndices(nextAssembly, nextReversed, currentFullSeqLength, nextAssemblyRefBases.length());

            currentFullSeqLength += nextAssemblyRefBases.length();

            if(overlapLength > 0)
            {
                sequenceCigar.append(format("%dI", overlapLength));
                mSequenceOverlaps.put(currentFullSeqLength - 1, link.overlapBases());
            }

            sequenceCigar.append(format("%dM", nextAssembly.refBaseLength()));

            if(hasOuterExtensions && i == assemblyCount - 1)
            {
                // add on the extension sequence for the last assembly
                String assemblyExtensionBases = nextAssembly.formJunctionSequence();

                fullSequence.append(assemblyExtensionBases);

                setAssemblyReadIndices(assembly, assemblyReversed, currentFullSeqLength, 0);

                currentFullSeqLength = assemblyExtensionBases.length();

                sequenceCigar.append(format("%dS", assemblyExtensionBases.length()));
            }

            lastAddedReversed = nextReversed;
        }

        mSequenceCigar = sequenceCigar.toString();

        return fullSequence.toString();
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

        for(SupportRead support : assembly.support())
        {
            int junctionReadStartDistance = support.junctionReadStartDistance();

            Orientation readOrientation = isReversed ? Orientation.REVERSE : Orientation.FORWARD;

            int fullSequenceIndex = existingSequenceLength;

            if(includeAssemblyRefBaseLength)
                fullSequenceIndex += assemblyRefBaseLength;

            if(!isReversed)
            {
                // typical scenario - eg for a junction read, this will bring the read's start index to earlier (lower) than the assembly index
                fullSequenceIndex -= junctionReadStartDistance;
            }
            else
            {
                // eg a junction read with JRSD = 100, ie the read starting 100 bases before the junction, read end at +50 into extension bases
                // after reversing, the read starts (from its end at -50 relative to this assembly's reversed junction end, extending to +100 past it
                fullSequenceIndex += junctionReadStartDistance - support.baseLength() + 1;
            }

            support.setFullAssemblyInfo(fullSequenceIndex, readOrientation);
        }

        // TODO: any secondaries need their support reads' full sequence assembly indices to be set too

    }

    public String assemblyCigar() { return mSequenceCigar; }

    public String toString()
    {
        return format("%s length(%d) breakends(%d)", info(), mFullSequenceLength, mBreakends.size());
    }

    // TODO: clean-up
    // public JunctionAssembly first() { return mAssemblies.get(0); }
    // public JunctionAssembly second() { return mAssemblies.size() > 1 ? mAssemblies.get(1) : null; }

    /*
    public StructuralVariantType svType()
    {
        return mAssemblies.size() == 1 ? StructuralVariantType.SGL : mAssemblyLink.svType();
    }

    public int svLength()
    {
        return mAssemblies.size() == 1 ? 0 : mAssemblyLink.length();
    }

    public int refBaseLength()
    {
        int refBaseLength = mAssemblies.get(0).refBaseLength();

        if(mAssemblies.size() > 1)
            refBaseLength += mAssemblies.get(1).refBaseLength();

        return refBaseLength;
    }
   */

    /* OLD FULL SEQUENCE BUILDING

        // factor in orientation and overlapping or inserted bases
        JunctionAssembly first, second;
        boolean firstReversed = false;
        boolean secondReversed = false;

        if(mAssemblyLink.first().junction().Orient != mAssemblyLink.second().junction().Orient)
        {
            if(mAssemblyLink.first().junction().isForward())
            {
                first = mAssemblyLink.first();
                second = mAssemblyLink.second();
            }
            else
            {
                first = mAssemblyLink.second();
                second = mAssemblyLink.first();
            }
        }
        else
        {
            first = mAssemblyLink.first();
            second = mAssemblyLink.second();

            if(mAssemblyLink.first().junction().isForward())
                secondReversed = true;
            else
                firstReversed = true;
        }

        // now build a sequence of first followed by second
        int overlapLength = mAssemblyLink.overlapBases().length();
        String insertedBases = mAssemblyLink.insertedBases();

        int fullSequenceLength = first.refBaseLength() + second.refBaseLength() + insertedBases.length() - overlapLength;

        StringBuilder fullSequence = new StringBuilder(fullSequenceLength);

        fullSequence.append(firstReversed ? Nucleotides.reverseComplementBases(first.formRefBaseSequence()) : first.formRefBaseSequence());

        if(!insertedBases.isEmpty())
        {
            if(firstReversed)
                fullSequence.append(Nucleotides.reverseComplementBases(insertedBases)); // keep inserted bases in the same direction as a full sequence
            else
                fullSequence.append(insertedBases);
        }

        String secondSequence = secondReversed ? Nucleotides.reverseComplementBases(second.formRefBaseSequence()) : second.formRefBaseSequence();

        if(overlapLength >= secondSequence.length())
            return null;

        fullSequence.append(overlapLength > 0 ? secondSequence.substring(overlapLength) : secondSequence);

        // set full assembly indices for all supporting reads, so they can be tested after alignment
        for(SupportRead read : first.support())
        {
            if(!firstReversed)
            {
                read.setLinkedAssemblyIndex(read.junctionAssemblyIndex());
            }
            else
            {
                int readEndIndex = read.junctionAssemblyIndex() + read.baseLength() - 1;
                int invertedStartIndex = first.baseLength() - readEndIndex - 1;
                read.setLinkedAssemblyIndex(invertedStartIndex);
            }
        }

        int secondStartAdjustment = first.refBaseLength() + insertedBases.length() - overlapLength;

        for(SupportRead read : second.support())
        {
            if(!secondReversed)
            {
                // use the read's relative position to its assembly's junction to set its position in the full assembly
                read.setLinkedAssemblyIndex(secondStartAdjustment - read.junctionReadIndex());
            }
            else
            {
                // it was index 3 vs length 10 (0-9), then inverted is 10 - 3 - 1 = 6
                int invertedReadJunctionIndex = read.baseLength() - read.junctionAssemblyIndex() - 1;
                read.setLinkedAssemblyIndex(secondStartAdjustment - invertedReadJunctionIndex);
            }
        }

        return fullSequence.toString();
        */
}
