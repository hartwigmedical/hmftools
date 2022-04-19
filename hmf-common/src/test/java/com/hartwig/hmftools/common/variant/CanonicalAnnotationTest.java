package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactoryTest;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.snpeff.ImmutableSnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CanonicalAnnotationTest {

    private final DriverGenePanel genePanel = DriverGenePanelFactoryTest.testGenePanel();
    private final Set<String> driverGenes = genePanel.driverGenes().stream().map(DriverGene::gene).collect(Collectors.toSet());
    private final List<HmfTranscriptRegion> transcripts = HmfGenePanelSupplier.allGeneList37();

    @Test
    public void testTrimEnsembleTranscriptId() {
        assertEquals("ENST00000361570", CanonicalAnnotation.trimEnsembleVersion("ENST00000361570"));
        assertEquals("ENST00000361570", CanonicalAnnotation.trimEnsembleVersion("ENST00000361570.v8"));
    }
}
