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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class HlaAlignment
{

    private final BwaMemAlignment baseAlignment;
    public int Position;
    private final int MapQuality;
    @NotNull private final String Cigar;

    public static @NotNull Set<HlaAlignment> hlaAlignments(BwaMemAlignment alignment, RefGenomeVersion refGenomeVersion)
    {
        Set<HlaAlignment> result = new HashSet<>();
        result.add(new HlaAlignment(alignment)); // TODO only if hla
        if (alignment.getMapQual() > 30)
        {
            return result;
        }
        if (alignment.getXATag() != null) {
            List<AlternativeAlignment> alternatives = AlternativeAlignment.fromLocationTag(alignment.getXATag());
            alternatives.stream()
                    .filter(a -> isHla(a,refGenomeVersion ))
                    .forEach(a -> result.add(new HlaAlignment(alignment, a)));
        }
        if (result.isEmpty())
        {
            throw new IllegalStateException("No HLA alignments found");
        }
        return result;
    }

    public HlaAlignment(final BwaMemAlignment baseAlignment, AlternativeAlignment alignment)
    {
        this.baseAlignment = baseAlignment;
        // TODO check in chr6
        Position = alignment.Position;
        // TODO  check in HLA
        MapQuality = 0; // alignment.MapQual; No. AlternativeAlignment.MapQual seems actually to be the edit distance
        Cigar = alignment.Cigar;
    }

    public HlaAlignment(final BwaMemAlignment baseAlignment)
    {
        this.baseAlignment = baseAlignment;
        // TODO check in chr6
        Position = baseAlignment.getRefStart() + 1;
        // TODO  check in HLA
        MapQuality = baseAlignment.getMapQual();
        Cigar = baseAlignment.getCigar();
    }

    public SAMRecord createSamRecord(SAMFileHeader header, RawFastaData raw, HlaAlignment mate)
    {
        SAMRecord remappedRecord = new SAMRecord(header);
        remappedRecord.setReadName(raw.ReadName);
        remappedRecord.setHeader(header);
        remappedRecord.setFlags(getSamFlag());
        remappedRecord.setReferenceIndex(getRefId());
        remappedRecord.setAlignmentStart(getRefStart());
        remappedRecord.setMappingQuality(getMapQual());
        remappedRecord.setCigarString(getCigar());
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
        remappedRecord.setAttribute(SamRecordUtils.NUM_MUTATONS_ATTRIBUTE, getNMismatches());
        remappedRecord.setAttribute(SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE, getMDTag());
        remappedRecord.setAttribute(SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE, getAlignerScore());
        remappedRecord.setAttribute(SamRecordUtils.SUBOPTIMAL_SCORE_ATTRIBUTE, getSuboptimalScore());

        remappedRecord.setMateReferenceIndex(mate.getRefId());
        remappedRecord.setMateAlignmentStart(mate.getRefStart());
        remappedRecord.setInferredInsertSize(calculateInsertSize(remappedRecord, mate));
        remappedRecord.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, mate.Cigar);
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
        return Position + TextCigarCodec.decode(Cigar).getReferenceLength() - 1;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final HlaAlignment that = (HlaAlignment) o;
        return Position == that.Position && Objects.equals(baseAlignment, that.baseAlignment);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(baseAlignment, Position);
    }

    @Override
    public String toString()
    {
        return "HlaAlignment{" +
                ", Position=" + Position +
                ", MapQuality=" + MapQuality +
                ", Cigar='" + Cigar + '\'' +
                '}';
    }

    public int getSamFlag()
    {
        return baseAlignment.getSamFlag();
    }

    public int getRefId()
    {
        return baseAlignment.getRefId();
    }

    public int getRefStart()
    {
        return Position;
    }

    public int getMapQual()
    {
        return MapQuality;
    }

    @NotNull
    public String getCigar()
    {
        return Cigar;
    }

    public Object getNMismatches()
    {
        return baseAlignment.getNMismatches();
    }

    public Object getMDTag()
    {
        return baseAlignment.getMDTag();
    }

    public Object getAlignerScore()
    {
        return baseAlignment.getAlignerScore();
    }

    public Object getSuboptimalScore()
    {
        return baseAlignment.getSuboptimalScore();
    }

    private static boolean isHla(AlternativeAlignment alternativeAlignment, RefGenomeVersion refGenomeVersion)
    {
        final List<ChrBaseRegion> hlaRegions = ImmuneRegions.getHlaRegions(refGenomeVersion);
        return hlaRegions.stream().anyMatch(chrBaseRegion -> chrBaseRegion.containsPosition(alternativeAlignment.Position));
    }
}
