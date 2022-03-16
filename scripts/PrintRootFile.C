
{
  std::string sipFile = "2018_SingleMu/L1T_JetMET_Res_def_nVtxgt50_part0_of_1.root";
  TFile *f=new TFile(sipFile.c_str());
  if ( ! f->IsOpen() )
  {
    printf("%s couldn't open\n",sipFile.c_str());
    return;
  }

  printf("File %s:\n",sipFile.c_str());
  //std::cout << f->ls() << std::endl;
  f->ls();

}
