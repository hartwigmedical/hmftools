package com.hartwig.hmftools.civic.data;

import java.util.List;

public abstract class CivicEvidenceItem {
    public abstract String id();

    public abstract String name();

    public abstract String description();

    public abstract String diseaseId();

    public abstract String diseaseName();

    public abstract String level();

    public abstract String significance();

    public abstract String direction();

    public abstract String origin();

    public abstract List<String> drugs();

}
