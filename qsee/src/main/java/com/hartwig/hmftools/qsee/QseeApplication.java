package com.hartwig.hmftools.qsee;

import static com.hartwig.hmftools.qsee.common.QseeConstants.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.qsee.plot.QseePlot;
import com.hartwig.hmftools.qsee.prep.QseePrep;
import com.hartwig.hmftools.qsee.prep.QseePrepConfig;

public class QseeApplication
{
    private final QseePrepConfig mConfig;

    public QseeApplication(QseePrepConfig config)
    {
        mConfig = config;
    }

    public void run()
    {
        new QseePrep(mConfig).run();
        new QseePlot(mConfig).run();
    }

    public static void main(String[] args){
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        QseePrepConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        QseePrepConfig config = new QseePrepConfig(configBuilder);
        new QseeApplication(config).run();
    }
}
