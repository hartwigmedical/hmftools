
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - SV Tools"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.svtools.simulation.ShatteringSim")
}
