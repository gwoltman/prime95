SECTIONS
{
  .data : { . = ALIGN(0x80); *(_GWDATA) }
}
