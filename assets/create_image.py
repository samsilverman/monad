from __future__ import annotations
from pathlib import Path
from PIL import Image


def create_side_by_side_image(input1: str, input2: str, output: str) -> None:
    image1 = Image.open(fp=input1).convert(mode='RGBA')
    image2 = Image.open(fp=input2).convert(mode='RGBA')

    # Combine horizontally
    width = 2 * image1.width
    height = image1.height

    combined_image = Image.new(mode='RGBA', size=(width, height))

    combined_image.paste(im=image1, box=(0, 0))
    combined_image.paste(im=image2, box=(image1.width, 0))

    combined_image.save(fp=output, format="PNG")


def create_triple_image(input1: str, input2: str, input3: str, output: str) -> None:
    image1 = Image.open(fp=input1).convert(mode='RGBA')
    image2 = Image.open(fp=input2).convert(mode='RGBA')
    image3 = Image.open(fp=input3).convert(mode='RGBA')

    # Combine horizontally (2x2 grid)
    width = 2 * image1.width
    height = 2 * image1.height

    combined_image = Image.new(mode='RGBA', size=(width, height))

    combined_image.paste(im=image1, box=(0, image1.height // 2))
    combined_image.paste(im=image2, box=(image1.width, 0))
    combined_image.paste(im=image3, box=(image1.width, image1.height))

    combined_image.save(fp=output, format="PNG")


def main() -> None:
    folder = Path(__file__).parent.resolve()

    input1 = folder / 'images'/ '4_LinearElasticity_1.png'
    input2 = folder / 'images'/ '4_LinearElasticity_2.png'
    output = folder / 'images'/ '4_LinearElasticity.png'

    create_side_by_side_image(input1=input1, input2=input2, output=output)

    input1 = folder / 'images'/ '5_LinearDielectric_1.png'
    input2 = folder / 'images'/ '5_LinearDielectric_2.png'
    output = folder / 'images'/ '5_LinearDielectric.png'

    create_side_by_side_image(input1=input1, input2=input2, output=output)

    input1 = folder / 'images'/ '6_LinearPiezoelectricity_1.png'
    input2 = folder / 'images'/ '6_LinearPiezoelectricity_2.png'
    input3 = folder / 'images'/ '6_LinearPiezoelectricity_3.png'
    output = folder / 'images'/ '6_LinearPiezoelectricity.png'

    create_triple_image(input1=input1, input2=input2, input3=input3, output=output)


if __name__ == '__main__':
    main()
