package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.Driver;

public interface GenomeWide extends Driver {
    enum Type {
        MSI,
        HRD,
        TMB_HIGH,
        TMB_LOW,
    }
}
