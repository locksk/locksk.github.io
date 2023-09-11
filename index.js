// main js stuff, i think relates to header responsiveness - part of initial template
$(document).ready(function(){
  $("a").on('click', function(event) {
    if (this.hash !== "") {
      event.preventDefault();
      var hash = this.hash;
      $('body,html').animate({
      scrollTop: $(hash).offset().top
      }, 1200, function(){
      window.location.hash = hash;
     });
     } 
    });
});

var width = $(window).width(); 

window.onscroll = function(){
if ((width >= 900)){
    if(document.body.scrollTop > 80 || document.documentElement.scrollTop > 80) {
        $("#middle").css("background-size","150% auto");
    }else{
        $("#middle").css("background-size","100% auto");        
    }
}
};

setTimeout(function(){
    $("#loading").addClass("animated fadeOut");
    setTimeout(function(){
      $("#loading").removeClass("animated fadeOut");
      $("#loading").css("display","none");
    },800);
},1450);

// scroll to top button 
var scrollToTopBtn = document.querySelector(".scrollToTopBtn");
var rootElement = document.documentElement;

function handleScroll() {
  // Do something on scroll
  var scrollTotal = rootElement.scrollHeight - rootElement.clientHeight;
  if (rootElement.scrollTop / scrollTotal > 0.18) {
    // Show button
    scrollToTopBtn.classList.add("showBtn");
  } else {
    // Hide button
    scrollToTopBtn.classList.remove("showBtn");
  }
}

function scrollToTop() {
  // Scroll to top logic
  rootElement.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}
scrollToTopBtn.addEventListener("click", scrollToTop);
document.addEventListener("scroll", handleScroll);

// tool tip stay

//$('.tooltip .bottom').on('mouseover',function(){
//  $(this).css({
//      "visibility": "visible",
//      "opacity": "1"
//     // more styles go here
//  });
//});

// lightbox code 
const images = [
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/DSC01759.jpg?raw=true', title: 'Aonach Eagach' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/DSC01257.jpg?raw=true', title: 'Schiehallion' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/DSC01896.jpg?raw=true', title: 'Three Sisters of Glen Coe' },
  // Add more images here
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/934A7860.jpg?raw=true', title: 'Oli Sykes of BMTH' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/IMG_0726.jpg?raw=true', title: 'James Arthur' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/IMG_0393.jpg?raw=true', title: 'Steve Morse of Deep Purple' },
  // Add more images here
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/Untitled_Panorama1 small.jpg?raw=true', title: 'Sgwd Gwladys' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/e934A0241.jpg?raw=true', title: 'Dinefwr Castle' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/DSC09575s.jpg?raw=true', title: 'Plunge pool at Sgwd Gwladys' },
  // Add more images here
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/DSC02441.jpg?raw=true', title: 'Model 1 & 2 Wall Clocks at Paulin Watches' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/IMG_8056.jpg?raw=true', title: 'Liam Gallagher, 2017' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/vertical_pano.jpg?raw=true', title: 'Milky way rising over Dunraven Bay' },
  // Add more images here
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/IMG_9044.jpg?raw=true', title: 'Hayley Williams of Paramore' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/eIMG_6183.jpg?raw=true', title: 'Penderyn Distillery' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/DSC02404.jpg?raw=true', title: 'Finishing Touches at Paulin Watches' },
  // Add more images here
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/e934A2235.jpg?raw=true', title: 'Lighthouse at Nash Point' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/eIMG_7966.jpg?raw=true', title: 'Anthony Joshua, 2017' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/IMG_8787 s.jpg?raw=true', title:'Le Palais du Justice en Nantes' },
  // Add more images here
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/ee934A3395.jpg?raw=true', title: 'Stars over Marloes Sands' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/IMG_2343.jpg?raw=true', title: 'Falls of Bruar' },
  { src: 'https://github.com/locksk/landing-page/blob/main/photos/DSC05786.jpg?raw=true', title: 'Shaping at Waterbourne Surfboards' },
  // Add more images here
];

let currentIndex = 0;

function openLightbox(imageSrc, title) {
    currentIndex = images.findIndex((image) => image.src === imageSrc);
    showImage(currentIndex);
    document.getElementById('lightbox').style.display = 'flex';

    // Add event listener to close on outside click
    document.addEventListener('click', outsideClick);
    // Add event listener to close on Esc key press
    document.addEventListener('keydown', escKeyPress);
    // Add event listener to navigate with left and right arrow keys
    document.addEventListener('keydown', navigateWithArrows);
}

function closeLightbox() {
    document.getElementById('lightbox').style.display = 'none';

    // Remove event listeners
    document.removeEventListener('click', outsideClick);
    document.removeEventListener('keydown', escKeyPress);
    document.removeEventListener('keydown', navigateWithArrows);
}

function outsideClick(event) {
    if (event.target === document.getElementById('lightbox')) {
        closeLightbox();
    }
}

function escKeyPress(event) {
    if (event.key === 'Escape') {
        closeLightbox();
    }
}

function navigateWithArrows(event) {
    if (event.key === 'ArrowLeft') {
        prevImage();
    } else if (event.key === 'ArrowRight') {
        nextImage();
    }
}

function showImage(index) {
    const image = images[index];
    document.getElementById('lightbox-image').src = image.src;
    document.getElementById('lightbox-title').textContent = image.title;
}

function prevImage() {
    if (currentIndex > 0) {
        currentIndex--;
        showImage(currentIndex);
    }
}

function nextImage() {
    if (currentIndex < images.length - 1) {
        currentIndex++;
        showImage(currentIndex);
    }
}