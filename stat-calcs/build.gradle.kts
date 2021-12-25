
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - Statistical Calcs"
version = "1.0"

dependencies {
    implementation(project(":hmf-common"))
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.statcalcs.cooc.CoOccurenceCalcs")
}
