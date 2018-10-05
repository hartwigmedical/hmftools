package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneModel {

    @NotNull
    public abstract List<HmfTranscriptRegion> somaticVariantGenePanel();
    @NotNull
    public abstract List<HmfTranscriptRegion> cnvGenePanel();
    @NotNull
    public abstract Map<String, DriverCategory> geneDriverCategoryMap();
    @NotNull
    public abstract Set<String> drupActionableGenes();

    @Value.Derived
    @Nullable
    public DriverCategory geneDriverCategory(@NotNull String gene) {
        return geneDriverCategoryMap().get(gene);
    }

    @Value.Derived
    public long somaticVariantsNumberOfBases() {
        return somaticVariantGenePanel().stream().mapToLong(GenomeRegion::bases).sum();
    }

    @Value.Derived
    public int somaticVariantNumberOfRegions() {
        return somaticVariantGenePanel().size();
    }
}
