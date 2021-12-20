package com.hartwig.hmftools.serve.refgenome;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RefGenomeResource {

    @NotNull
    public abstract IndexedFastaSequenceFile refSequence();

    @NotNull
    public abstract List<DriverGene> driverGenes();

    @NotNull
    public abstract KnownFusionCache knownFusionCache();

    @NotNull
    public abstract EnsemblDataCache ensemblDataCache();

    @NotNull
    public abstract Map<RefGenomeVersion, String> chainToOtherRefGenomeMap();

    @NotNull
    public abstract ProteinResolver proteinResolver();

}
