package com.hartwig.hmftools.qsee;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.plot.QseePlot;
import com.hartwig.hmftools.qsee.plot.QseePlotConfig;
import com.hartwig.hmftools.qsee.prep.QseePrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

public class QseeApplication
{
    ConfigBuilder mConfig;

    public QseeApplication(ConfigBuilder configBuilder)
    {
        mConfig = configBuilder;
    }

    public void run()
    {
        QseePrepConfig prepConfig = new QseePrepConfig(mConfig);
        QseePrep qseePrep = new QseePrep(prepConfig);
        qseePrep.run();

        QseePlotConfig plotConfig = new QseePlotConfig(mConfig);
        QseePlot qseePlot = new QseePlot(plotConfig);
        qseePlot.run();
    }

    public static void main(String[] args){
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseePrepConfig.registerConfig(configBuilder);
        QseePlotConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        new QseeApplication(configBuilder).run();
    }
}
