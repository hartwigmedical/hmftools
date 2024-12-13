package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.immune.ImmuneRegions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

public class PreferredAlignment
{
    @NotNull private final BwaMemAlignment baseAlignment;
    @NotNull private final RefGenomeVersion refGenomeVersion;
    private int RefStart;
    private int MapQuality;
    @NotNull private String Cigar;

    public PreferredAlignment(@NotNull final BwaMemAlignment baseAlignment, final BwaMemAlignment mateAlignment, @NotNull final RefGenomeVersion refGenomeVersion)
    {
        this.baseAlignment = baseAlignment;
        this.refGenomeVersion = refGenomeVersion;
        final List<PossibleAlignment> hlaAlternatives = alternativeAlignments()
                .stream()
                .filter(this::isHla)
                .filter(alternativeAlignment -> PossibleAlignment.areDirectionallyCompatible(mateAlignment, alternativeAlignment))
                .map(alternativeAlignment -> new PossibleAlignment(mateAlignment, alternativeAlignment))
                .sorted().collect(Collectors.toList());
        if(!hlaAlternatives.isEmpty())
        {
            final PossibleAlignment bestAlternative = hlaAlternatives.get(0);
            // If the base alignment is on chr6 then it needs to be included as
            // an alternative as it may be closer to the mate than any alternative.
            if (baseAlignment.getRefId() == 5)
            {
                String chromosome = bestAlternative.alignment().Chromosome;
                AlternativeAlignment baseAlternative = new AlternativeAlignment(chromosome, baseAlignment.getRefStart(), bestAlternative.alignment().Orient, baseAlignment.getCigar(), baseAlignment.getMapQual());
                PossibleAlignment basePossibleAlignment = new PossibleAlignment(mateAlignment, baseAlternative);
                if (basePossibleAlignment.distanceFromMate() <= bestAlternative.distanceFromMate())
                {
                    setFromBase();
                } else {
                    RefStart = bestAlternative.alignment().Position;
                    MapQuality = bestAlternative.alignment().MapQual;
                    Cigar = bestAlternative.alignment().Cigar;
                }
            }
        }
        else
        {
            setFromBase();
        }
    }

    private void setFromBase()
    {
        RefStart = baseAlignment.getRefStart() + 1;
        MapQuality = baseAlignment.getMapQual();
        Cigar = baseAlignment.getCigar();
    }

    private List<AlternativeAlignment> alternativeAlignments()
    {
        return baseAlignment.getXATag() == null ? List.of() : AlternativeAlignment.fromLocationTag(baseAlignment.getXATag());
    }

    private boolean isHla(AlternativeAlignment alternativeAlignment)
    {
        final List<ChrBaseRegion> hlaRegions = ImmuneRegions.getHlaRegions(refGenomeVersion);
        return hlaRegions.stream().anyMatch(chrBaseRegion -> chrBaseRegion.containsPosition(alternativeAlignment.Position));
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
        return RefStart;
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
}
