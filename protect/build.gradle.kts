
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - PROTECT"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":serve"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
    testImplementation(project(path = ":hmf-common", configuration = "tests"))
    testImplementation(project(path = ":serve", configuration = "tests"))
}

application {
    mainClass.set("com.hartwig.hmftools.protect.ProtectApplication")
}
