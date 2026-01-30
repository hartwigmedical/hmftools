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

import org.umccr.java.hellbender.utils.bwa.BwaMemAlignment;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class HlaAlignment
{

    private final BwaMemAlignment BaseAlignment;
    public final int Position;
    final Set<SAMFlag> Flags;
    private final int mMapQuality;
    private final String mCigar;

    public static Set<HlaAlignment> hlaAlignments(BwaMemAlignment alignment, RefGenomeVersion refGenomeVersion)
    {
        Set<HlaAlignment> result = new HashSet<>();
        result.add(new HlaAlignment(alignment));
        if(alignment.getMapQual() > 30)
        {
            return result;
        }
        if(alignment.getXATag() != null)
        {
            List<AltAlignment> alternatives = AltAlignment.fromLocationTag(alignment.getXATag());
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

    public HlaAlignment(final BwaMemAlignment baseAlignment, AltAlignment alignment)
    {
        BaseAlignment = baseAlignment;
        Position = alignment.Position;
        mMapQuality = 0; // Only the edit distance is available in the alts in the xa tag.
        mCigar = alignment.Cigar;
        Flags = SAMFlag.getFlags(getSamFlag());
    }

    public HlaAlignment(final BwaMemAlignment baseAlignment)
    {
        this.BaseAlignment = baseAlignment;
        Position = baseAlignment.getRefStart() + 1;
        mMapQuality = baseAlignment.getMapQual();
        mCigar = baseAlignment.getCigar();
        Flags = SAMFlag.getFlags(getSamFlag());
    }

    public SAMRecord createSamRecord(SAMFileHeader header, BamReadData raw, HlaAlignment mate)
    {
        SAMRecord remappedRecord = new SAMRecord(header);
        remappedRecord.setReadName(raw.ReadName);
        remappedRecord.setHeader(header);
        remappedRecord.setFlags(getSamFlag());
        if(isUnmapped())
        {
            remappedRecord.setReferenceIndex(mate.getRefId());
            remappedRecord.setAlignmentStart(mate.Position);
            remappedRecord.setCigarString("*");
            remappedRecord.setMappingQuality(0);
        }
        else
        {
            remappedRecord.setReferenceIndex(getRefId());
            remappedRecord.setAlignmentStart(Position);
            remappedRecord.setCigarString(mCigar);
            remappedRecord.setMappingQuality(mMapQuality);
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
        remappedRecord.setAttribute(SamRecordUtils.XS_ATTRIBUTE, BaseAlignment.getSuboptimalScore());

        if(mate.isUnmapped())
        {
            remappedRecord.setMateAlignmentStart(Position);
            remappedRecord.setMateReferenceIndex(getRefId());
            remappedRecord.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, "*");
        }
        else
        {
            remappedRecord.setMateAlignmentStart(mate.Position);
            remappedRecord.setMateReferenceIndex(mate.getRefId());
            remappedRecord.setAttribute(SamRecordUtils.MATE_CIGAR_ATTRIBUTE, mate.mCigar);
        }
        remappedRecord.setInferredInsertSize(calculateInsertSize(remappedRecord, mate));
        remappedRecord.setAttribute(SamRecordUtils.MATE_QUALITY_ATTRIBUTE, mate.mMapQuality);

        return remappedRecord;
    }

    private static int calculateInsertSize(final SAMRecord record, final HlaAlignment mate)
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
        return Position + TextCigarCodec.decode(mCigar).getReferenceLength() - 1;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final HlaAlignment that = (HlaAlignment) o;
        return Position == that.Position && Objects.equals(BaseAlignment, that.BaseAlignment);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(BaseAlignment, Position);
    }

    @Override
    public String toString()
    {
        return "HlaAlignment{" +
                ", Position=" + Position +
                ", MapQuality=" + mMapQuality +
                ", Cigar='" + mCigar + '\'' +
                '}';
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

    private static boolean isHla(AltAlignment alternativeAlignment, RefGenomeVersion refGenomeVersion)
    {
        final List<ChrBaseRegion> hlaRegions = ImmuneRegions.getHlaRegions(refGenomeVersion);
        return hlaRegions.stream().anyMatch(chrBaseRegion -> chrBaseRegion.containsPosition(alternativeAlignment.Position));
    }
}
