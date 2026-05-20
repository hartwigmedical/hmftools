package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.ValidationStringency;

public class CompareUtils
{
    @VisibleForTesting
    public static List<String> compareReads(final SAMRecord orig, final SAMRecord newRead, final CompareConfig config)
    {
        if(orig.getReadUnmappedFlag() && newRead.getReadUnmappedFlag())
            return Collections.emptyList();

        boolean unmappedDiff = orig.getReadUnmappedFlag() != newRead.getReadUnmappedFlag();
        boolean mateUnmappedDiff = orig.getMateUnmappedFlag() != newRead.getMateUnmappedFlag();
        boolean skipUnmapping = config.IgnoreReduxUnmapped && (unmappedDiff || mateUnmappedDiff);

        boolean ownCigarsEquiv = config.CigarBoundaryTolerance > 0
                && !orig.getReadUnmappedFlag() && !newRead.getReadUnmappedFlag()
                && cigarsEquivalent(orig, newRead, config.CigarBoundaryTolerance);

        boolean mateCigarsEquiv = config.CigarBoundaryTolerance > 0
                && mateCigarsEquivalent(orig, newRead, config.CigarBoundaryTolerance);

        List<String> diffs = new ArrayList<>();

        if(!skipUnmapping && !ownCigarsEquiv && coordsDiffer(orig, newRead))
            diffs.add(format("coords(%s/%s)", readCoords(orig), readCoords(newRead)));

        if(config.CompareCoordsOnly)
            return diffs;

        if(!skipUnmapping && !ownCigarsEquiv && !mateCigarsEquiv
                && orig.getInferredInsertSize() != newRead.getInferredInsertSize())
            diffs.add(format("insertSize(%d/%d)", orig.getInferredInsertSize(), newRead.getInferredInsertSize()));

        if(!skipUnmapping && !ownCigarsEquiv && !orig.getCigarString().equals(newRead.getCigarString()))
            diffs.add(format("cigar(%s/%s)", orig.getCigarString(), newRead.getCigarString()));

        if(!skipUnmapping && orig.getMappingQuality() != newRead.getMappingQuality())
            diffs.add(format("mapQuality(%d/%d)", orig.getMappingQuality(), newRead.getMappingQuality()));

        if(orig.getFlags() != newRead.getFlags())
            addFlagDiffs(orig, newRead, diffs, config, skipUnmapping, unmappedDiff, mateUnmappedDiff);

        if(!skipUnmapping && orig.getReadPairedFlag() && newRead.getReadPairedFlag()
                && !orig.getMateUnmappedFlag() && !newRead.getMateUnmappedFlag())
            addMateFieldDiffs(orig, newRead, diffs, config.CigarBoundaryTolerance);

        addAttributeDiffs(orig, newRead, diffs, config, skipUnmapping, mateCigarsEquiv);

        if(config.CheckBasesAndQuals)
            addBasesAndQualsDiffs(orig, newRead, diffs);

        return diffs;
    }

    private static boolean coordsDiffer(final SAMRecord orig, final SAMRecord newRead)
    {
        if(orig.getAlignmentStart() != newRead.getAlignmentStart())
            return true;
        if(orig.getAlignmentEnd() != newRead.getAlignmentEnd())
            return true;
        return orig.getContig() != null && newRead.getContig() != null
                && !orig.getContig().equals(newRead.getContig());
    }

    private static void addFlagDiffs(final SAMRecord orig, final SAMRecord newRead, final List<String> diffs,
            final CompareConfig config, final boolean skipUnmapping, final boolean unmappedDiff, final boolean mateUnmappedDiff)
    {
        // strand flag is undefined on an unmapped record per SAM spec
        if(!skipUnmapping && !unmappedDiff && orig.getReadNegativeStrandFlag() != newRead.getReadNegativeStrandFlag())
            diffs.add(format("negStrand(%s/%s)", orig.getReadNegativeStrandFlag(), newRead.getReadNegativeStrandFlag()));

        if(!config.IgnoreDupDiffs && orig.getDuplicateReadFlag() != newRead.getDuplicateReadFlag())
            diffs.add(format("duplicate(%s/%s)", orig.getDuplicateReadFlag(), newRead.getDuplicateReadFlag()));

        if(skipUnmapping)
            return;

        if(unmappedDiff)
            diffs.add(format("unmapped(%s/%s)", orig.getReadUnmappedFlag(), newRead.getReadUnmappedFlag()));

        if(mateUnmappedDiff)
            diffs.add(format("mateUnmapped(%s/%s)", orig.getMateUnmappedFlag(), newRead.getMateUnmappedFlag()));
    }

    private static void addMateFieldDiffs(final SAMRecord orig, final SAMRecord newRead, final List<String> diffs, final int tolerance)
    {
        if(!Objects.equals(orig.getMateReferenceName(), newRead.getMateReferenceName()))
            diffs.add(format("mateChr(%s/%s)", orig.getMateReferenceName(), newRead.getMateReferenceName()));

        int pnextDelta = Math.abs(orig.getMateAlignmentStart() - newRead.getMateAlignmentStart());
        if(pnextDelta > 0 && pnextDelta > tolerance)
            diffs.add(format("matePos(%d/%d)", orig.getMateAlignmentStart(), newRead.getMateAlignmentStart()));

        if(orig.getMateNegativeStrandFlag() != newRead.getMateNegativeStrandFlag())
            diffs.add(format("mateNegStrand(%s/%s)", orig.getMateNegativeStrandFlag(), newRead.getMateNegativeStrandFlag()));
    }

    private static void addAttributeDiffs(final SAMRecord orig, final SAMRecord newRead, final List<String> diffs,
            final CompareConfig config, final boolean skipUnmapping, final boolean mateCigarsEquiv)
    {
        if(!config.IgnoreSupplementaryAttribute
                && !orig.hasAttribute(UNMAP_ATTRIBUTE) && !newRead.hasAttribute(UNMAP_ATTRIBUTE))
            addAttribDiff(orig, newRead, diffs, SUPPLEMENTARY_ATTRIBUTE);

        if(!skipUnmapping && !mateCigarsEquiv)
            addAttribDiff(orig, newRead, diffs, MATE_CIGAR_ATTRIBUTE);
    }

    private static void addAttribDiff(final SAMRecord orig, final SAMRecord newRead, final List<String> diffs, final String attr)
    {
        String origValue = orig.getStringAttribute(attr);
        String newValue = newRead.getStringAttribute(attr);
        if(!Objects.equals(origValue, newValue))
            diffs.add(format("attrib_%s(%s/%s)", attr,
                    origValue == null ? "missing" : origValue,
                    newValue == null ? "missing" : newValue));
    }

    private static void addBasesAndQualsDiffs(final SAMRecord orig, final SAMRecord newRead, final List<String> diffs)
    {
        if(!orig.getReadString().equals(newRead.getReadString()))
            diffs.add(format("bases(%s/%s)", orig.getReadString(), newRead.getReadString()));

        if(!orig.getBaseQualityString().equals(newRead.getBaseQualityString()))
            diffs.add(format("baseQual(%s/%s)", orig.getBaseQualityString(), newRead.getBaseQualityString()));
    }

    private static String readCoords(final SAMRecord read)
    {
        if(!read.getReadUnmappedFlag())
            return format("%s:%d-%d:%s", read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(), read.getCigarString());
        return format("%s:%d", read.getContig(), read.getAlignmentStart());
    }

    static boolean cigarsEquivalent(final SAMRecord orig, final SAMRecord newRead, final int toleranceBp)
    {
        if(orig.getReadNegativeStrandFlag() != newRead.getReadNegativeStrandFlag())
            return false;
        if(orig.getContig() == null || newRead.getContig() == null || !orig.getContig().equals(newRead.getContig()))
            return false;
        // alignments must share a start or end within tolerance — otherwise the cigars belong to two
        // distinct placements that happen to have similar element shapes
        if(Math.abs(orig.getAlignmentStart() - newRead.getAlignmentStart()) > toleranceBp
                && Math.abs(orig.getAlignmentEnd() - newRead.getAlignmentEnd()) > toleranceBp)
            return false;
        return cigarElementsEquivalent(orig.getCigar(), newRead.getCigar(), toleranceBp);
    }

    static boolean mateCigarsEquivalent(final SAMRecord orig, final SAMRecord newRead, final int toleranceBp)
    {
        if(!orig.getReadPairedFlag() || !newRead.getReadPairedFlag())
            return false;
        if(orig.getMateUnmappedFlag() || newRead.getMateUnmappedFlag())
            return false;
        if(Math.abs(orig.getMateAlignmentStart() - newRead.getMateAlignmentStart()) > toleranceBp)
            return false;
        if(orig.getMateNegativeStrandFlag() != newRead.getMateNegativeStrandFlag())
            return false;
        if(!Objects.equals(orig.getMateReferenceName(), newRead.getMateReferenceName()))
            return false;

        String origMc = orig.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        String newMc = newRead.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        if(origMc == null || newMc == null)
            return false;
        if(origMc.equals(newMc))
            return true;
        return cigarElementsEquivalent(TextCigarCodec.decode(origMc), TextCigarCodec.decode(newMc), toleranceBp);
    }

    private static boolean cigarElementsEquivalent(final Cigar a, final Cigar b, final int toleranceBp)
    {
        List<CigarElement> elementsA = absorbEndSoftClips(a.getCigarElements(), toleranceBp);
        List<CigarElement> elementsB = absorbEndSoftClips(b.getCigarElements(), toleranceBp);
        if(elementsA.size() != elementsB.size() || elementsA.isEmpty())
            return false;

        for(int i = 0; i < elementsA.size(); ++i)
        {
            CigarElement ea = elementsA.get(i);
            CigarElement eb = elementsB.get(i);
            if(ea.getOperator() != eb.getOperator())
                return false;
            if(Math.abs(ea.getLength() - eb.getLength()) > toleranceBp)
                return false;
        }
        return true;
    }

    // merge a leading or trailing soft-clip of length <= toleranceBp into its adjacent M-class op
    private static List<CigarElement> absorbEndSoftClips(final List<CigarElement> elements, final int toleranceBp)
    {
        if(elements.size() < 2)
            return elements;

        List<CigarElement> result = new ArrayList<>(elements);
        int last = result.size() - 1;
        if(isAbsorbableSoftClip(result.get(last), toleranceBp) && isMatchClass(result.get(last - 1).getOperator()))
        {
            CigarElement neighbour = result.get(last - 1);
            result.set(last - 1, new CigarElement(neighbour.getLength() + result.get(last).getLength(), neighbour.getOperator()));
            result.remove(last);
        }

        if(result.size() >= 2 && isAbsorbableSoftClip(result.get(0), toleranceBp) && isMatchClass(result.get(1).getOperator()))
        {
            CigarElement neighbour = result.get(1);
            result.set(1, new CigarElement(neighbour.getLength() + result.get(0).getLength(), neighbour.getOperator()));
            result.remove(0);
        }
        return result;
    }

    private static boolean isAbsorbableSoftClip(final CigarElement element, final int toleranceBp)
    {
        return element.getOperator() == CigarOperator.S && element.getLength() <= toleranceBp;
    }

    private static boolean isMatchClass(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X;
    }

    public static SamReaderFactory makeSamReaderFactory(final CompareConfig config)
    {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT);
        if(config.RefGenomeFile != null && !config.RefGenomeFile.isEmpty())
            readerFactory.referenceSequence(new File(config.RefGenomeFile));
        return readerFactory;
    }
}
