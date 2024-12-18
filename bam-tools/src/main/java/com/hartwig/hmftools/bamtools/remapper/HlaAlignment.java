package com.hartwig.hmftools.bamtools.remapper;

import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.immune.ImmuneRegions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

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
        if (alignment.getXATag() != null) {
            List<AlternativeAlignment> alternatives = AlternativeAlignment.fromLocationTag(alignment.getXATag());
            alternatives.stream()
                    .filter(a -> isHla(a,refGenomeVersion ))
                    .forEach(a -> {result.add(new HlaAlignment(alignment, a));});
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
        MapQuality = alignment.MapQual;
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
