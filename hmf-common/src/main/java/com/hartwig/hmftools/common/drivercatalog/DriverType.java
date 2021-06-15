package com.hartwig.hmftools.common.drivercatalog;

import java.util.EnumSet;

public enum DriverType
{
    AMP,
    PARTIAL_AMP,
    DEL,
    HOM_DISRUPTION,
    MUTATION,
    GERMLINE;

    public static final EnumSet<DriverType> DRIVERS_GERMLINE = EnumSet.of(DriverType.GERMLINE);
    public static final EnumSet<DriverType> DRIVERS_LINX = EnumSet.of(DriverType.HOM_DISRUPTION);
    public static final EnumSet<DriverType> DRIVERS_COPY_NUMBER = EnumSet.of(DriverType.AMP, DriverType.PARTIAL_AMP, DriverType.DEL);
    public static final EnumSet<DriverType> DRIVERS_MUTATION = EnumSet.of(DriverType.MUTATION);

}
