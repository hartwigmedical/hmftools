package com.hartwig.hmftools.serve.extraction;

import java.util.List;

import com.hartwig.hmftools.common.serve.datamodel.characteristic.TumorCharacteristic;
import com.hartwig.hmftools.common.serve.datamodel.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.immuno.ImmunoHLA;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EventExtractorOutput {

    @Nullable
    public abstract List<VariantHotspot> hotspots();

    @Nullable
    public abstract List<CodonAnnotation> codons();

    @Nullable
    public abstract List<ExonAnnotation> exons();

    @Nullable
    public abstract GeneLevelAnnotation geneLevelEvent();

    @Nullable
    public abstract KnownCopyNumber knownCopyNumber();

    @Nullable
    public abstract KnownFusionPair knownFusionPair();

    @Nullable
    public abstract TumorCharacteristic characteristic();

    @Nullable
    public abstract ImmunoHLA hla();
}
