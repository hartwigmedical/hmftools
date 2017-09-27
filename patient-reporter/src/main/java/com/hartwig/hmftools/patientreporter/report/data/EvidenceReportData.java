package com.hartwig.hmftools.patientreporter.report.data;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.jasperreports.engine.data.JRBeanCollectionDataSource;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class EvidenceReportData {
    public abstract List<Alteration> getAlterations();

    @Value.Derived
    public List<Alteration> getAlterationsWithEvidence() {
        return getAlterations().stream().filter(alteration -> alteration.getEvidence().size() > 0).collect(Collectors.toList());
    }

    @NotNull
    public JRBeanCollectionDataSource toDataSource() {
        return new JRBeanCollectionDataSource(Lists.newArrayList(this));
    }
}
