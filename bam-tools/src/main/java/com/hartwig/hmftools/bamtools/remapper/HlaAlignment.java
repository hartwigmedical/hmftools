package com.hartwig.hmftools.bamtools.remapper;

import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.immune.ImmuneRegions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Arrays;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class HlaAlignment
{

    private final BwaMemAlignment BaseAlignment;
    public final int Position_1Based;
    private final int MapQuality;
    @NotNull
    private final String Cigar;
    @NotNull
    final Set<SAMFlag> Flags;

    public static @NotNull Set<HlaAlignment> hlaAlignments(BwaMemAlignment alignment, RefGenomeVersion refGenomeVersion)
    {
        Set<HlaAlignment> result = new HashSet<>();
        result.add(new HlaAlignment(alignment)); // TODO only if hla
        if(alignment.getMapQual() > 30)
        {
            return result;
        }
        if(alignment.getXATag() != null)
        {
            List<AlternativeAlignment> alternatives = AlternativeAlignment.fromLocationTag(alignment.getXATag());
            alternatives.stream()
                    .filter(a -> isHla(a, refGenomeVersion))//
                    .forEach(a -> result.add(new HlaAlignment(alignment, a)));
        }
        if(result.isEmpty())
        {
            throw new IllegalStateException("No HLA alignments found");
        }
        return result;
    }

    public HlaAlignment(final BwaMemAlignment baseAlignment, AlternativeAlignment alignment)
    {
        this.BaseAlignment = baseAlignment;
        // TODO check in chr6
        Position_1Based = alignment.Position;
        // TODO  check in HLA
        MapQuality = 0; // alignment.MapQual; No. AlternativeAlignment.MapQual seems actually to be the edit distance
        Cigar = alignment.Cigar;
        Flags = SAMFlag.getFlags(getSamFlag());
    }

    public HlaAlignment(final BwaMemAlignment baseAlignment)
    {
        this.BaseAlignment = baseAlignment;
        // TODO check in chr6
        Position_1Based = baseAlignment.getRefStart() + 1;
        // TODO  check in HLA
        MapQuality = baseAlignment.getMapQual();
        Cigar = baseAlignment.getCigar();
        Flags = SAMFlag.getFlags(getSamFlag());
    }

    public SAMRecord createSamRecord(SAMFileHeader header, RawFastaData raw, HlaAlignment mate)
    {
        SAMRecord remappedRecord = new SAMRecord(header);
        remappedRecord.setReadName(raw.ReadName);
        remappedRecord.setHeader(header);
        remappedRecord.setFlags(getSamFlag());
        if(isUnmapped())
        {
            remappedRecord.setReferenceIndex(mate.getRefId());
            remappedRecord.setAlignmentStart(mate.Position_1Based);
            remappedRecord.setCigarString("*");
            remappedRecord.setMappingQuality(0);
        }
        else
        {
            remappedRecord.setReferenceIndex(getRefId());
            remappedRecord.setAlignmentStart(Position_1Based);
            remappedRecord.setCigarString(Cigar);
            remappedRecord.setMappingQuality(MapQuality);
        }
        if(remappedRecord.getReadNegativeStrandFlag())
        {
            remappedRecord.setReadBases(Nucleotides.reverseComplementBases(raw.Bases));
            remappedRecord.setBaseQualities(Arrays.reverseArray(raw.Qualities));
        }
        else
        {
            remappedRecord.setReadBases(raw.Bases);
            remappedRecord.setBaseQualities(raw.Qualities);
        }
        remappedRecord.setAttribute(SamRecordUtils.NUM_MUTATONS_ATTRIBUTE, BaseAlignment.getNMismatches());
        remappedRecord.setAttribute(SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE, BaseAlignment.getMDTag());
        remappedRecord.setAttribute(SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE, BaseAlignment.getAlignerScore());
        remappedRecord.setAttribute(SamRecordUtils.SUBOPTIMAL_SCORE_ATTRIBUTE, BaseAlignment.getSuboptimalScore());

        if (mate.isUnmapped())
        {
            remappedRecord.setMateAlignmentStart(Position_1Based);
            remappedRecord.setMateReferenceIndex(getRefId());
            remappedRecord.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, "*");
        }
        else
        {
            remappedRecord.setMateAlignmentStart(mate.Position_1Based);
            remappedRecord.setMateReferenceIndex(mate.getRefId());
            remappedRecord.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, mate.Cigar);
        }
        remappedRecord.setInferredInsertSize(calculateInsertSize(remappedRecord, mate));
        remappedRecord.setAttribute(SamRecordUtils.MATE_QUALITY_ATTRIBUTE, mate.MapQuality);

        return remappedRecord;
    }

    private static int calculateInsertSize(@NotNull final SAMRecord record, @NotNull final HlaAlignment mate)
    {
        if(record.isSecondaryOrSupplementary())
        {
            return 0;
        }
        if(!Objects.equals(record.getReferenceIndex(), mate.getRefId()))
        {
            return 0;
        }

        if(record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
        {
            return 0;
        }
        if(record.getReadNegativeStrandFlag())
        {
            return -1 * (record.getAlignmentEnd() - record.getMateAlignmentStart() + 1);
        }
        return mate.getAlignmentEnd() - record.getAlignmentStart() + 1;
    }

    private int getAlignmentEnd()
    {
        return Position_1Based + TextCigarCodec.decode(Cigar).getReferenceLength() - 1;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final HlaAlignment that = (HlaAlignment) o;
        return Position_1Based == that.Position_1Based && Objects.equals(BaseAlignment, that.BaseAlignment);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(BaseAlignment, Position_1Based);
    }

    @Override
    public String toString()
    {
        return "HlaAlignment{" +
                ", Position=" + Position_1Based +
                ", MapQuality=" + MapQuality +
                ", Cigar='" + Cigar + '\'' +
                '}';
    }

    public boolean isReadReverseStrand()
    {
        return Flags.contains(SAMFlag.READ_REVERSE_STRAND);
    }

    public boolean isUnmapped()
    {
        return Flags.contains(SAMFlag.READ_UNMAPPED);
    }

    public int getSamFlag()
    {
        return BaseAlignment.getSamFlag();
    }

    public Integer getRefId()
    {
        return BaseAlignment.getRefId();
    }

    private static boolean isHla(AlternativeAlignment alternativeAlignment, RefGenomeVersion refGenomeVersion)
    {
        final List<ChrBaseRegion> hlaRegions = ImmuneRegions.getHlaRegions(refGenomeVersion);
        return hlaRegions.stream().anyMatch(chrBaseRegion -> chrBaseRegion.containsPosition(alternativeAlignment.Position));
    }
}
