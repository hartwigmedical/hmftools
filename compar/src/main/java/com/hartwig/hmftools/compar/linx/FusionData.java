package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.Category.LINX_DATA;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.diffValue;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvAnnotation;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;

import org.apache.commons.compress.utils.Lists;

public class FusionData implements ComparableItem
{
    public final LinxFusion Fusion;

    public FusionData(final LinxFusion fusion)
    {
        Fusion = fusion;
    }

    public Category category() { return FUSION; }

    public boolean matches(final ComparableItem other)
    {
        final FusionData otherFusion = (FusionData)other;

        return otherFusion.Fusion.name().equals(Fusion.name());
    }

    public List<String> findDifferences(final ComparableItem other)
    {
        final FusionData otherFusion = (FusionData)other;

        final List<String> diffs = Lists.newArrayList();

        if(Fusion.reported() != otherFusion.Fusion.reported())
        {
            diffs.add(String.format("reported(%s/%s)", Fusion.reported(), otherFusion.Fusion.reported()));
        }

        if(!Fusion.reportedType().equals(otherFusion.Fusion.reportedType()))
        {
            diffs.add(String.format("reportedType(%s/%s)", Fusion.reportedType(), otherFusion.Fusion.reportedType()));
        }

        if(Fusion.phased() != otherFusion.Fusion.phased())
        {
            diffs.add(String.format("phased(%s/%s)", Fusion.phased(), otherFusion.Fusion.phased()));
        }

        if(Fusion.likelihood() != otherFusion.Fusion.likelihood())
        {
            diffs.add(String.format("likelihood(%s/%s)", Fusion.likelihood(), otherFusion.Fusion.likelihood()));
        }

        if(Fusion.fusedExonUp() != otherFusion.Fusion.fusedExonUp() || Fusion.fusedExonDown() != otherFusion.Fusion.fusedExonDown())
        {
            diffs.add(String.format("fusedExons(%d_%d/%d_%d)",
                    Fusion.fusedExonUp(), Fusion.fusedExonDown(), otherFusion.Fusion.fusedExonUp(), otherFusion.Fusion.fusedExonDown()));
        }

        /*
            public abstract int chainLength();
            public abstract int chainLinks();
            public abstract boolean chainTerminated();
            public abstract String domainsKept();
            public abstract String domainsLost();
            public abstract int skippedExonsUp();
            public abstract int skippedExonsDown();

         */

        return diffs;
    }

    public String description()
    {
        return String.format("%s", Fusion.name());
    }
}
