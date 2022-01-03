
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - Statistical Calcs"
version = "1.0"

dependencies {
    implementation(project(":hmf-common"))
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.statcalcs.cooc.CoOccurenceCalcs")
}
