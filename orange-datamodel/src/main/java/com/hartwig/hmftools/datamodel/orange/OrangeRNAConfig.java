package com.hartwig.hmftools.datamodel.orange;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRNAConfig {

    Logger LOGGER = LogManager.getLogger(OrangeRNAConfig.class);

    @NotNull
    String rnaSampleId();

    @NotNull
    String isofoxGeneDistributionCsv();

    @NotNull
    String isofoxAltSjCohortCsv();

    @NotNull
    String isofoxSummaryCsv();

    @NotNull
    String isofoxGeneDataCsv();

    @NotNull
    String isofoxFusionCsv();

    @NotNull
    String isofoxAltSpliceJunctionCsv();

}
